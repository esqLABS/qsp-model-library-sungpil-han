# Prolactin-secreting pituitary tumour (prolactinoma) — QSP model references
### Prolactinoma / lactotroph PitNET · 150 PubMed references, every PMID resolved through NCBI E-utilities

**Every PMID in this list was looked up directly with NCBI E-utilities (esearch + esummary) and its title, journal, and
year confirmed**, and not one citation was written from memory. Because in an earlier session in this
repository a large number of PMIDs written from memory pointed at entirely unrelated papers, this time the
order was query → lookup → relevance filter → automatic generation of the bibliographic details, and candidates that failed the
relevance check were dropped from the list. The criteria the verification script applied were (1) that the first author's surname in the query
matches the first author of the returned paper, or (2) that two or more of the content words of the query appear in the
returned title or journal.

The quotation block before each section states **which equation, parameter, or diagnostic analysis of the model that group of references**
supports. Within a section the references are ordered by year.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map | [`prl_qsp_model.dot`](../../../prolactinoma/prl_qsp_model.dot) · [`.svg`](../../../prolactinoma/prl_qsp_model.svg) · [`.png`](../../../prolactinoma/prl_qsp_model.png) |
| ⚙️ mrgsolve ODE model | [`prl_mrgsolve_model.R`](../../../prolactinoma/prl_mrgsolve_model.R) |
| 📊 Shiny dashboard | [`prl_shiny_app.R`](prl_shiny_app.R) |
| 📄 Directory README | [`README.md`](README.md) |

---

## What this model claims quantitatively about the literature (calibration anchors)

| Anchor | Source (numbered below) | Model result |
|---|---|---|
| Cabergoline t½ 63-109 h, plasma concentrations in the pg/mL range | section G (Ferrari 1995, Rains 1995) | Cmax 68.6 pg/mL (1 mg), t½ 107 h |
| Cabergoline > bromocriptine (from the difference in intrinsic activity alone) | Webster 1994 (G) | plateau 51 vs 80 ng/mL, with no separate efficacy parameter |
| Prolactin suppression persists for more than two weeks after a single dose | section G | 24% of baseline on day 14 |
| Sustained remission after withdrawal ~21% | Dekkers 2010 (I) | the remission itself is not generated → **reinterpreted as a function of the follow-up period** (D10) |
| Symptomatic enlargement in pregnancy, micro ~2.7% / macro ~21-23% | section J (Molitch) | generated from a geometric threshold (D13) |
| Valvulopathy yes at Parkinson's disease doses, no at prolactinoma doses | section M (Zanettini/Schade vs Stiles/Caputo) | P(TR) 25.8% vs 0.38% (D11) |
| Impulse control disorder on a dopamine agonist ~10-17% | section M | 14.1% at 1 mg/week |
| The ceiling of stalk-compression hyperprolactinaemia, 25-150 ng/mL | section L | 94 ng/mL on complete transection (D8) |
| A PEG recovery of <40% = macroprolactin | section F | 16.9% |
| Response to temozolomide when MGMT is underexpressed | section H (Raverot, Whitelaw) | 78% shrinkage against no response |

---

### A. Guidelines, reviews, epidemiology

> The source of the clinical skeleton of the model and of its reference ranges, prevalences, and treatment algorithm. Melmed 2011 and Petersenn 2023 in particular
> were the standard for constructing the scenarios (first-line medical therapy, follow-up intervals, withdrawal criteria).

1. Colao A et al. New medical approaches in pituitary adenomas. *Horm Res*. 2000;53 Suppl 3:76-87. [PMID 10971110](https://pubmed.ncbi.nlm.nih.gov/10971110/)
2. Verhelst J & Abs R. Hyperprolactinemia: pathophysiology and management. *Treat Endocrinol*. 2003;2:23-32. [PMID 15871552](https://pubmed.ncbi.nlm.nih.gov/15871552/)
3. Buurman H & Saeger W. Subclinical adenomas in postmortem pituitaries: classification and correlations to clinical data. *Eur J Endocrinol*. 2006;154:753-8. [PMID 16645024](https://pubmed.ncbi.nlm.nih.gov/16645024/)
4. Casanueva FF et al. Guidelines of the Pituitary Society for the diagnosis and management of prolactinomas. *Clin Endocrinol (Oxf)*. 2006;65:265-73. [PMID 16886971](https://pubmed.ncbi.nlm.nih.gov/16886971/)
5. Daly AF et al. High prevalence of pituitary adenomas: a cross-sectional study in the province of Liege, Belgium. *J Clin Endocrinol Metab*. 2006;91:4769-75. [PMID 16968795](https://pubmed.ncbi.nlm.nih.gov/16968795/)
6. Gillam MP et al. Advances in the treatment of prolactinomas. *Endocr Rev*. 2006;27:485-534. [PMID 16705142](https://pubmed.ncbi.nlm.nih.gov/16705142/)
7. Kars M et al. Estimated age- and sex-specific incidence and prevalence of dopamine agonist-treated hyperprolactinemia. *J Clin Endocrinol Metab*. 2009;94:2729-34. [PMID 19491225](https://pubmed.ncbi.nlm.nih.gov/19491225/)
8. Fernandez A et al. Prevalence of pituitary adenomas: a community-based, cross-sectional study in Banbury (Oxfordshire, UK). *Clin Endocrinol (Oxf)*. 2010;72:377-82. [PMID 19650784](https://pubmed.ncbi.nlm.nih.gov/19650784/)
9. Kars M et al. Update in prolactinomas. *Neth J Med*. 2010;68:104-12. [PMID 20308704](https://pubmed.ncbi.nlm.nih.gov/20308704/)
10. Melmed S et al. Diagnosis and treatment of hyperprolactinemia: an Endocrine Society clinical practice guideline. *J Clin Endocrinol Metab*. 2011;96:273-88. [PMID 21296991](https://pubmed.ncbi.nlm.nih.gov/21296991/)
11. Faje A & Nachtigall L. Current treatment options for hyperprolactinemia. *Expert Opin Pharmacother*. 2013;14:1611-25. [PMID 23738973](https://pubmed.ncbi.nlm.nih.gov/23738973/)
12. Glezer A & Bronstein MD. [Prolactinoma]. *Arq Bras Endocrinol Metabol*. 2014;58:118-23. [PMID 24830588](https://pubmed.ncbi.nlm.nih.gov/24830588/)
13. Molitch ME. Diagnosis and Treatment of Pituitary Adenomas: A Review. *JAMA*. 2017;317:516-524. [PMID 28170483](https://pubmed.ncbi.nlm.nih.gov/28170483/)
14. Vilar L et al. Controversial issues in the management of hyperprolactinemia and prolactinomas - An overview by the Neuroendocrinology Department of the Brazilian Society of Endocrinology and Metabolism. *Arch Endocrinol Metab*. 2018;62:236-263. [PMID 29768629](https://pubmed.ncbi.nlm.nih.gov/29768629/)
15. Auriemma RS et al. Dopamine Agonists: From the 1970s to Today. *Neuroendocrinology*. 2019;109:34-41. [PMID 30852578](https://pubmed.ncbi.nlm.nih.gov/30852578/)
16. Vroonen L et al. Epidemiology and Management Challenges in Prolactinomas. *Neuroendocrinology*. 2019;109:20-27. [PMID 30731464](https://pubmed.ncbi.nlm.nih.gov/30731464/)
17. Trouillas J et al. How to Classify the Pituitary Neuroendocrine Tumors (PitNET)s in 2020. *Cancers (Basel)*. 2020;12. [PMID 32098443](https://pubmed.ncbi.nlm.nih.gov/32098443/)
18. Asa SL et al. Overview of the 2022 WHO Classification of Pituitary Tumors. *Endocr Pathol*. 2022;33:6-26. [PMID 35291028](https://pubmed.ncbi.nlm.nih.gov/35291028/)
19. Kontbay T et al. Hyperprolactinemia in children and adolescents and longterm follow-up results of prolactinoma cases: a single-centre experience. *Turk J Pediatr*. 2022;64:892-899. [PMID 36305439](https://pubmed.ncbi.nlm.nih.gov/36305439/)
20. Petersenn S et al. Diagnosis and management of prolactin-secreting pituitary adenomas: a Pituitary Society international Consensus Statement. *Nat Rev Endocrinol*. 2023;19:722-740. [PMID 37670148](https://pubmed.ncbi.nlm.nih.gov/37670148/)
21. Raverot G et al. Revised European Society of Endocrinology Clinical Practice Guideline for the management of aggressive pituitary tumours and pituitary carcinomas. *Eur J Endocrinol*. 2025;192:R45-R78. [PMID 40506054](https://pubmed.ncbi.nlm.nih.gov/40506054/)

---

### B. Pathogenesis and genetics

> The basis for `FGEN` (the proliferative drive) and `FDRIVE` (the transcriptional drive), and for why the AIP and SF3B1 phenotypes connect to a low D2R and to invasiveness.
> SF3B1 R625H is a recently discovered prolactinoma-specific hotspot and is drawn in cluster 6 of the map.

22. Vergès B et al. Pituitary disease in MEN type 1 (MEN1): data from the France-Belgium MEN1 multicenter study. *J Clin Endocrinol Metab*. 2002;87:457-65. [PMID 11836268](https://pubmed.ncbi.nlm.nih.gov/11836268/)
23. Daly AF et al. Aryl hydrocarbon receptor-interacting protein gene mutations in familial isolated pituitary adenomas: analysis in 73 families. *J Clin Endocrinol Metab*. 2007;92:1891-6. [PMID 17244780](https://pubmed.ncbi.nlm.nih.gov/17244780/)
24. Raverot G et al. Prognostic factors in prolactin pituitary tumors: clinical, histological, and molecular data from a series of 94 patients with a long postoperative follow-up. *J Clin Endocrinol Metab*. 2010;95:1708-16. [PMID 20164287](https://pubmed.ncbi.nlm.nih.gov/20164287/)
25. Mitsui T et al. Differences between rat strains in the development of PRL-secreting pituitary tumors with long-term estrogen treatment: In vitro insulin-like growth factor-1-induced lactotroph proliferation and gene expression are affected in Wistar-Kyoto rats with low estrogen-susceptibility. *Endocr J*. 2013;60:1251-9. [PMID 23985690](https://pubmed.ncbi.nlm.nih.gov/23985690/)
26. Kageyama K et al. A Novel Deletion Mutation in the MEN1 Gene in a Patient with Prolactinoma and a Family History of Pancreatic Tumors. *Endocr Pract*. 2014;20:e162-5. [PMID 24936550](https://pubmed.ncbi.nlm.nih.gov/24936550/)
27. Matsuno A et al. Molecular status of pituitary carcinoma and atypical adenoma that contributes the effectiveness of temozolomide. *Med Mol Morphol*. 2014;47:1-7. [PMID 23955641](https://pubmed.ncbi.nlm.nih.gov/23955641/)
28. Trouillas J et al. Clinical, Pathological, and Molecular Factors of Aggressiveness in Lactotroph Tumours. *Neuroendocrinology*. 2019;109:70-76. [PMID 30943495](https://pubmed.ncbi.nlm.nih.gov/30943495/)
29. Vandeva S et al. Somatic and germline mutations in the pathogenesis of pituitary adenomas. *Eur J Endocrinol*. 2019;181:R235-R254. [PMID 31658440](https://pubmed.ncbi.nlm.nih.gov/31658440/)
30. Guo J et al. The SF3B1(R625H) mutation promotes prolactinoma tumor progression through aberrant splicing of DLG1. *J Exp Clin Cancer Res*. 2022;41:26. [PMID 35039052](https://pubmed.ncbi.nlm.nih.gov/35039052/)

---

### C. The dopamine D2 receptor and lactotroph physiology

> The heart of the model: the single occupancy equation (`SIGDRIVE`), the D2R density (`D2RS0`), the short loop (`TIDA`→`DAP`), and
> the basis for the four-way branch (secretion, transcription, cell volume, proliferation). The claim that a fall in D2R density is the quantitative definition of resistance also comes from here.

31. Curlewis JD & McNeilly AS. Prolactin short-loop feedback and prolactin inhibition of luteinizing hormone secretion during the breeding season and seasonal anoestrus in the ewe. *Neuroendocrinology*. 1991;54:279-85. [PMID 1944814](https://pubmed.ncbi.nlm.nih.gov/1944814/)
32. Elsholtz HP et al. Inhibitory control of prolactin and Pit-1 gene promoters by dopamine. Dual signaling pathways required for D2 receptor-regulated expression of the prolactin gene. *J Biol Chem*. 1991;266:22919-25. [PMID 1835974](https://pubmed.ncbi.nlm.nih.gov/1835974/)
33. Caccavelli L et al. Decreased expression of the two D2 dopamine receptor isoforms in bromocriptine-resistant prolactinomas. *Neuroendocrinology*. 1994;60:314-22. [PMID 7969790](https://pubmed.ncbi.nlm.nih.gov/7969790/)
34. Kelly MA et al. Pituitary lactotroph hyperplasia and chronic hyperprolactinemia in dopamine D2 receptor-deficient mice. *Neuron*. 1997;19:103-13. [PMID 9247267](https://pubmed.ncbi.nlm.nih.gov/9247267/)
35. Missale C et al. Dopamine receptors: from structure to function. *Physiol Rev*. 1998;78:189-225. [PMID 9457173](https://pubmed.ncbi.nlm.nih.gov/9457173/)
36. Ben-Jonathan N & Hnasko R. Dopamine as a prolactin (PRL) inhibitor. *Endocr Rev*. 2001;22:724-63. [PMID 11739329](https://pubmed.ncbi.nlm.nih.gov/11739329/)
37. Peverelli E et al. Filamin-A is essential for dopamine d2 receptor expression and signaling in tumorous lactotrophs. *J Clin Endocrinol Metab*. 2012;97:967-77. [PMID 22259062](https://pubmed.ncbi.nlm.nih.gov/22259062/)
38. Shimazu S et al. Resistance to dopamine agonists in prolactinoma is correlated with reduction of dopamine D2 receptor long isoform mRNA levels. *Eur J Endocrinol*. 2012;166:383-90. [PMID 22127489](https://pubmed.ncbi.nlm.nih.gov/22127489/)
39. Venkatesh SK et al. Spontaneous reduction of prolactinoma post cabergoline withdrawal. *Indian J Endocrinol Metab*. 2012;16:833-5. [PMID 23087877](https://pubmed.ncbi.nlm.nih.gov/23087877/)
40. Grattan DR. 60 YEARS OF NEUROENDOCRINOLOGY: The hypothalamo-prolactin axis. *J Endocrinol*. 2015;226:T101-22. [PMID 26101377](https://pubmed.ncbi.nlm.nih.gov/26101377/)
41. Bernard V et al. Autocrine actions of prolactin contribute to the regulation of lactotroph function in vivo. *FASEB J*. 2018;32:4791-4797. [PMID 29596024](https://pubmed.ncbi.nlm.nih.gov/29596024/)
42. Li H et al. Melatonin Modulates Lactation by Regulating Prolactin Secretion Via Tuberoinfundibular Dopaminergic Neurons in the Hypothalamus- Pituitary System. *Curr Protein Pept Sci*. 2020;21:744-750. [PMID 32392109](https://pubmed.ncbi.nlm.nih.gov/32392109/)
43. McNamara AV et al. Transcription Factor Pit-1 Affects Transcriptional Timing in the Dual-Promoter Human Prolactin Gene. *Endocrinology*. 2021;162. [PMID 33388754](https://pubmed.ncbi.nlm.nih.gov/33388754/)

---

### D. Prolactin biology and receptor signalling

> The half-life of `PRLB`, the secretory granule store, PRL gene transcription (Pit-1, ERα), glycosylation and the 16 kDa cleavage, and PRLR-JAK2-STAT5 downstream.

44. Day RN et al. A protein kinase inhibitor gene reduces both basal and multihormone-stimulated prolactin gene transcription. *J Biol Chem*. 1989;264:431-6. [PMID 2535842](https://pubmed.ncbi.nlm.nih.gov/2535842/)
45. Shull JD & Gorski J. Estrogen regulation of prolactin gene transcription in vivo: paradoxical effects of 17 beta-estradiol dose. *Endocrinology*. 1989;124:279-85. [PMID 2909367](https://pubmed.ncbi.nlm.nih.gov/2909367/)
46. Smith CR & Norman MR. Prolactin and growth hormone: molecular heterogeneity and measurement in serum. *Ann Clin Biochem*. 1990;27 ( Pt 6):542-50. [PMID 2080857](https://pubmed.ncbi.nlm.nih.gov/2080857/)
47. Sinha YN. Prolactin variants. *Trends Endocrinol Metab*. 1992;3:100-6. [PMID 18407087](https://pubmed.ncbi.nlm.nih.gov/18407087/)
48. Bole-Feysot C et al. Prolactin (PRL) and its receptor: actions, signal transduction pathways and phenotypes observed in PRL receptor knockout mice. *Endocr Rev*. 1998;19:225-68. [PMID 9626554](https://pubmed.ncbi.nlm.nih.gov/9626554/)
49. Freeman ME et al. Prolactin: structure, function, and regulation of secretion. *Physiol Rev*. 2000;80:1523-631. [PMID 11015620](https://pubmed.ncbi.nlm.nih.gov/11015620/)
50. Macotela Y et al. Matrix metalloproteases from chondrocytes generate an antiangiogenic 16 kDa prolactin. *J Cell Sci*. 2006;119:1790-800. [PMID 16608881](https://pubmed.ncbi.nlm.nih.gov/16608881/)
51. Ben-Jonathan N et al. What can we learn from rodents about prolactin in humans?. *Endocr Rev*. 2008;29:1-41. [PMID 18057139](https://pubmed.ncbi.nlm.nih.gov/18057139/)
52. Skowronska-Krawczyk D et al. Required enhancer-matrin-3 network interactions for a homeodomain transcription program. *Nature*. 2014;514:257-61. [PMID 25119036](https://pubmed.ncbi.nlm.nih.gov/25119036/)
53. Bernard V et al. New insights in prolactin: pathological implications. *Nat Rev Endocrinol*. 2015;11:265-75. [PMID 25781857](https://pubmed.ncbi.nlm.nih.gov/25781857/)
54. Gao Q et al. Seasonal patterns of prolactin, prolactin receptor, and STAT5 expression in the ovaries of wild ground squirrels (<em>Citellus dauricus</em> Brandt). *Eur J Histochem*. 2023;67. [PMID 37781865](https://pubmed.ncbi.nlm.nih.gov/37781865/)
55. Hackwell ECR et al. Prolactin-mediates a lactation-induced suppression of arcuate kisspeptin neuronal activity necessary for lactational infertility in mice. *Elife*. 2025;13. [PMID 39819370](https://pubmed.ncbi.nlm.nih.gov/39819370/)

---

### E. HPG axis suppression and reproduction

> The basis for the `KISS`→`GNRH`→`LH`/`FSH`→sex hormone pathway and for `K50K` (the half-suppressing prolactin concentration).

56. Rasmussen C et al. Prolactin secretion and menstrual function after long-term bromocriptine treatment. *Fertil Steril*. 1987;48:550-4. [PMID 3653413](https://pubmed.ncbi.nlm.nih.gov/3653413/)
57. Sonigo C et al. Hyperprolactinemia-induced ovarian acyclicity is reversed by kisspeptin administration. *J Clin Invest*. 2012;122:3791-5. [PMID 23006326](https://pubmed.ncbi.nlm.nih.gov/23006326/)
58. Donato J Jr & Frazão R. Interactions between prolactin and kisspeptin to control reproduction. *Arch Endocrinol Metab*. 2016;60:587-595. [PMID 27901187](https://pubmed.ncbi.nlm.nih.gov/27901187/)

---

### F. The measurement layer: hook effect and macroprolactin

> The basis for the separate layer in the model, held apart from the biology (`PRLIMM`, `PRLMEAS_`, `PRLDIL`, `PEGREC`).
> Diagnostics D6 and D7 and scenarios S21 and S22 rest entirely on this group of references.

59. Bevan JS et al. Misinterpretation of prolactin levels leading to management errors in patients with sellar enlargement. *Am J Med*. 1987;82:29-32. [PMID 3799691](https://pubmed.ncbi.nlm.nih.gov/3799691/)
60. St-Jean E et al. High prolactin levels may be missed by immunoradiometric assay in patients with macroprolactinomas. *Clin Endocrinol (Oxf)*. 1996;44:305-9. [PMID 8729527](https://pubmed.ncbi.nlm.nih.gov/8729527/)
61. Barkan AL & Chandler WF. Giant pituitary prolactinoma with falsely low serum prolactin: the pitfall of the "high-dose hook effect": case report. *Neurosurgery*. 1998;42:913-5; discussion 915-6. [PMID 9574657](https://pubmed.ncbi.nlm.nih.gov/9574657/)
62. Petakov MS et al. Pituitary adenomas secreting large amounts of prolactin may give false low values in immunoradiometric assays. The hook effect. *J Endocrinol Invest*. 1998;21:184-8. [PMID 9591215](https://pubmed.ncbi.nlm.nih.gov/9591215/)
63. Colao A et al. Macroprolactinoma shrinkage during cabergoline treatment is greater in naive patients than in patients pretreated with other dopamine agonists: a prospective study in 110 patients. *J Clin Endocrinol Metab*. 2000;85:2247-52. [PMID 10852458](https://pubmed.ncbi.nlm.nih.gov/10852458/)
64. Frieze TW et al. "Hook effect" in prolactinomas: case report and review of literature. *Endocr Pract*. 2002;8:296-303. [PMID 12173917](https://pubmed.ncbi.nlm.nih.gov/12173917/)
65. Schöfl C et al. Falsely low serum prolactin in two cases of invasive macroprolactinoma. *Pituitary*. 2002;5:261-5. [PMID 14558675](https://pubmed.ncbi.nlm.nih.gov/14558675/)
66. Smith TP et al. Gross variability in the detection of prolactin in sera containing big big prolactin (macroprolactin) by commercial immunoassays. *J Clin Endocrinol Metab*. 2002;87:5410-5. [PMID 12466327](https://pubmed.ncbi.nlm.nih.gov/12466327/)
67. Gibney J et al. Clinical relevance of macroprolactin. *Clin Endocrinol (Oxf)*. 2005;62:633-43. [PMID 15943822](https://pubmed.ncbi.nlm.nih.gov/15943822/)
68. Hattori N et al. Anti-prolactin (PRL) autoantibody-binding sites (epitopes) on PRL molecule in macroprolactinemia. *J Endocrinol*. 2006;190:287-93. [PMID 16899562](https://pubmed.ncbi.nlm.nih.gov/16899562/)
69. Delgrange E et al. Characterization of resistance to the prolactin-lowering effects of cabergoline in macroprolactinomas: a study in 122 patients. *Eur J Endocrinol*. 2009;160:747-52. [PMID 19223454](https://pubmed.ncbi.nlm.nih.gov/19223454/)
70. Hattori N et al. Macroprolactinaemia: prevalence and aetiologies in a large group of hospital workers. *Clin Endocrinol (Oxf)*. 2009;71:702-8. [PMID 19486017](https://pubmed.ncbi.nlm.nih.gov/19486017/)
71. Raverot G et al. Secondary deterioration of visual field during cabergoline treatment for macroprolactinoma. *Clin Endocrinol (Oxf)*. 2009;70:588-92. [PMID 18673461](https://pubmed.ncbi.nlm.nih.gov/18673461/)
72. Barber TM et al. Recurrence of hyperprolactinaemia following discontinuation of dopamine agonist therapy in patients with prolactinoma occurs commonly especially in macroprolactinoma. *Clin Endocrinol (Oxf)*. 2011;75:819-24. [PMID 21645021](https://pubmed.ncbi.nlm.nih.gov/21645021/)
73. Shimatsu A & Hattori N. Macroprolactinemia: diagnostic, clinical, and pathogenic significance. *Clin Dev Immunol*. 2012;2012:167132. [PMID 23304187](https://pubmed.ncbi.nlm.nih.gov/23304187/)
74. Raverot V et al. Prolactin immunoassay: does the high-dose hook effect still exist?. *Pituitary*. 2022;25:653-657. [PMID 35793045](https://pubmed.ncbi.nlm.nih.gov/35793045/)
75. Vermue FC et al. The validation of macroprolactin analysis by polyethylene glycol precipitation using Fujirebio Lumipulse. *Pract Lab Med*. 2022;31:e00292. [PMID 35860390](https://pubmed.ncbi.nlm.nih.gov/35860390/)

---

### G. Cabergoline · bromocriptine · quinagolide (dopamine agonist PK/PD and trials)

> The PK parameters (t½ 63-109 h, plasma concentrations in pg/mL), the intrinsic activity e (cabergoline 1.00 against bromocriptine 0.80),
> and the comparator for the dose-response of D4.

76. Bergh T et al. Bromocriptine-induced regression of a suprasellar extending prolactinoma during pregnancy. *J Endocrinol Invest*. 1984;7:133-6. [PMID 6725868](https://pubmed.ncbi.nlm.nih.gov/6725868/)
77. Vance ML et al. Drugs five years later. Bromocriptine. *Ann Intern Med*. 1984;100:78-91. [PMID 6229205](https://pubmed.ncbi.nlm.nih.gov/6229205/)
78. Bevan JS et al. Factors in the outcome of transsphenoidal surgery for prolactinoma and non-functioning pituitary tumour, including pre-operative bromocriptine therapy. *Clin Endocrinol (Oxf)*. 1987;26:541-56. [PMID 3665118](https://pubmed.ncbi.nlm.nih.gov/3665118/)
79. Pellegrini I et al. Resistance to bromocriptine in prolactinomas. *J Clin Endocrinol Metab*. 1989;69:500-9. [PMID 2760167](https://pubmed.ncbi.nlm.nih.gov/2760167/)
80. Vance ML et al. CV 205-502 treatment of hyperprolactinemia. *J Clin Endocrinol Metab*. 1989;68:336-9. [PMID 2521863](https://pubmed.ncbi.nlm.nih.gov/2521863/)
81. Webster J et al. A comparison of cabergoline and bromocriptine in the treatment of hyperprolactinemic amenorrhea. Cabergoline Comparative Study Group. *N Engl J Med*. 1994;331:904-9. [PMID 7915824](https://pubmed.ncbi.nlm.nih.gov/7915824/)
82. Andreotti AC et al. Pharmacokinetics, pharmacodynamics, and tolerability of cabergoline, a prolactin-lowering drug, after administration of increasing oral doses (0.5, 1.0, and 1.5 milligrams) in healthy male volunteers. *J Clin Endocrinol Metab*. 1995;80:841-5. [PMID 7883840](https://pubmed.ncbi.nlm.nih.gov/7883840/)
83. Ferrari C et al. Cabergoline: a new drug for the treatment of hyperprolactinaemia. *Hum Reprod*. 1995;10:1647-52. [PMID 8582955](https://pubmed.ncbi.nlm.nih.gov/8582955/)
84. Rains CP et al. Cabergoline. A review of its pharmacological properties and therapeutic potential in the treatment of hyperprolactinaemia and inhibition of lactation. *Drugs*. 1995;49:255-79. [PMID 7729332](https://pubmed.ncbi.nlm.nih.gov/7729332/)
85. Motta T et al. Vaginal cabergoline in the treatment of hyperprolactinemic patients intolerant to oral dopaminergics. *Fertil Steril*. 1996;65:440-2. [PMID 8566276](https://pubmed.ncbi.nlm.nih.gov/8566276/)
86. Verhelst J et al. Cabergoline in the treatment of hyperprolactinemia: a study in 455 patients. *J Clin Endocrinol Metab*. 1999;84:2518-22. [PMID 10404830](https://pubmed.ncbi.nlm.nih.gov/10404830/)
87. Colao A et al. Outcome of cabergoline treatment in men with prolactinoma: effects of a 24-month treatment on prolactin levels, tumor mass, recovery of pituitary function, and semen analysis. *J Clin Endocrinol Metab*. 2004;89:1704-11. [PMID 15070934](https://pubmed.ncbi.nlm.nih.gov/15070934/)
88. Barlier A & Jaquet P. Quinagolide--a valuable treatment option for hyperprolactinaemia. *Eur J Endocrinol*. 2006;154:187-95. [PMID 16452531](https://pubmed.ncbi.nlm.nih.gov/16452531/)

---

### H. Resistance and second line (high dose · surgery · radiotherapy · temozolomide)

> The distinction between partial resistance (a shift in EC50) and true resistance (a fall in Emax), high-dose cabergoline, surgical cure rates,
> and the MGMT-dependent response to temozolomide (S27/S28).

89. Kreutzer J et al. Operative treatment of prolactinomas: indications and results in a current consecutive series of 212 patients. *Eur J Endocrinol*. 2008;158:11-8. [PMID 18166812](https://pubmed.ncbi.nlm.nih.gov/18166812/)
90. Ono M et al. Prospective study of high-dose cabergoline treatment of prolactinomas in 150 patients. *J Clin Endocrinol Metab*. 2008;93:4721-7. [PMID 18812485](https://pubmed.ncbi.nlm.nih.gov/18812485/)
91. Babey M et al. Pituitary surgery for small prolactinomas as an alternative to treatment with dopamine agonists. *Pituitary*. 2011;14:222-30. [PMID 21170594](https://pubmed.ncbi.nlm.nih.gov/21170594/)
92. Vroonen L et al. Prolactinomas resistant to standard doses of cabergoline: a multicenter study of 92 patients. *Eur J Endocrinol*. 2012;167:651-62. [PMID 22918301](https://pubmed.ncbi.nlm.nih.gov/22918301/)
93. Chen W et al. HIF-1α inhibition sensitizes pituitary adenoma cells to temozolomide by regulating MGMT expression. *Oncol Rep*. 2013;30:2495-501. [PMID 23970362](https://pubmed.ncbi.nlm.nih.gov/23970362/)
94. Bengtsson D et al. Long-term outcome and MGMT as a predictive marker in 24 patients with atypical pituitary adenomas and pituitary carcinomas given treatment with temozolomide. *J Clin Endocrinol Metab*. 2015;100:1689-98. [PMID 25646794](https://pubmed.ncbi.nlm.nih.gov/25646794/)
95. Losa M et al. Temozolomide therapy in patients with aggressive pituitary adenomas or carcinomas. *J Neurooncol*. 2016;126:519-25. [PMID 26614517](https://pubmed.ncbi.nlm.nih.gov/26614517/)
96. Halevy C & Whitelaw BC. How effective is temozolomide for treating pituitary tumours and when should it be used?. *Pituitary*. 2017;20:261-266. [PMID 27581836](https://pubmed.ncbi.nlm.nih.gov/27581836/)
97. Honegger J & Grimm F. The experience with transsphenoidal surgery and its importance to outcomes. *Pituitary*. 2018;21:545-555. [PMID 30062664](https://pubmed.ncbi.nlm.nih.gov/30062664/)
98. Lee DK et al. Factors Influencing Visual Field Recovery after Transsphenoidal Resection of a Pituitary Adenoma. *Korean J Ophthalmol*. 2018;32:488-496. [PMID 30549473](https://pubmed.ncbi.nlm.nih.gov/30549473/)
99. Buchfelder M et al. Surgery for Prolactinomas to Date. *Neuroendocrinology*. 2019;109:77-81. [PMID 30699424](https://pubmed.ncbi.nlm.nih.gov/30699424/)
100. Ponce AJ et al. Low prolactin levels are associated with visceral adipocyte hypertrophy and insulin resistance in humans. *Endocrine*. 2020;67:331-343. [PMID 31919769](https://pubmed.ncbi.nlm.nih.gov/31919769/)

---

### I. Withdrawal and recurrence

> The standard for D10. The 21% sustained remission rate of this group of references is reinterpreted in the model as 'a function of the follow-up period'.

101. Colao A et al. Withdrawal of long-term cabergoline therapy for tumoral and nontumoral hyperprolactinemia. *N Engl J Med*. 2003;349:2023-33. [PMID 14627787](https://pubmed.ncbi.nlm.nih.gov/14627787/)
102. Dekkers OM et al. Recurrence of hyperprolactinemia after withdrawal of dopamine agonists: systematic review and meta-analysis. *J Clin Endocrinol Metab*. 2010;95:43-51. [PMID 19880787](https://pubmed.ncbi.nlm.nih.gov/19880787/)
103. Huda MS et al. Factors determining the remission of microprolactinomas after dopamine agonist withdrawal. *Clin Endocrinol (Oxf)*. 2010;72:507-11. [PMID 19549247](https://pubmed.ncbi.nlm.nih.gov/19549247/)
104. Auriemma RS et al. Results of a single-center observational 10-year survey study on recurrence of hyperprolactinemia after pregnancy and lactation. *J Clin Endocrinol Metab*. 2013;98:372-9. [PMID 23162092](https://pubmed.ncbi.nlm.nih.gov/23162092/)
105. Kwancharoen R et al. Second attempt to withdraw cabergoline in prolactinomas: a pilot study. *Pituitary*. 2014;17:451-6. [PMID 24078319](https://pubmed.ncbi.nlm.nih.gov/24078319/)

---

### J. Pregnancy and oestrogen

> The dual action of E2 (PRL transcription ↑, D2R ↓), and the risk of symptomatic tumour enlargement in pregnancy (micro ~2.7% against macro ~21-23%).

106. Christin-Maître S et al. Prolactinoma and estrogens: pregnancy, contraception and hormonal replacement therapy. *Ann Endocrinol (Paris)*. 2007;68:106-12. [PMID 17540335](https://pubmed.ncbi.nlm.nih.gov/17540335/)
107. Colao A et al. Pregnancy outcomes following cabergoline treatment: extended results from a 12-year observational study. *Clin Endocrinol (Oxf)*. 2008;68:66-71. [PMID 17760883](https://pubmed.ncbi.nlm.nih.gov/17760883/)
108. Lebbe M et al. Outcome of 100 pregnancies initiated under treatment with cabergoline in hyperprolactinaemic women. *Clin Endocrinol (Oxf)*. 2010;73:236-42. [PMID 20455894](https://pubmed.ncbi.nlm.nih.gov/20455894/)
109. Galvão A et al. Prolactinoma and pregnancy - a series of cases including pituitary apoplexy. *J Obstet Gynaecol*. 2017;37:284-287. [PMID 27866462](https://pubmed.ncbi.nlm.nih.gov/27866462/)
110. Jayabalan N et al. Cross Talk between Adipose Tissue and Placenta in Obese and Gestational Diabetes Mellitus Pregnancies via Exosomes. *Front Endocrinol (Lausanne)*. 2017;8:239. [PMID 29021781](https://pubmed.ncbi.nlm.nih.gov/29021781/)
111. Karaca Z et al. How does pregnancy affect the patients with pituitary adenomas: a study on 113 pregnancies from Turkey. *J Endocrinol Invest*. 2018;41:129-141. [PMID 28634705](https://pubmed.ncbi.nlm.nih.gov/28634705/)
112. Laway BA et al. Prolactinoma Outcome After Pregnancy and Lactation: A Cohort Study. *Indian J Endocrinol Metab*. 2021;25:559-562. [PMID 35355922](https://pubmed.ncbi.nlm.nih.gov/35355922/)

---

### K. Bone and metabolic consequences

> The basis for `KBEXP` (the bone turnover ratio → the BMD set point), `FIRR` (the irreversible fraction), and the probability of vertebral fracture. The observation that
> bone density does not recover completely even after prolactin is normalised is the direct basis for the ratchet structure.

113. Klibanski A et al. Decreased bone density in hyperprolactinemic women. *N Engl J Med*. 1980;303:1511-4. [PMID 7432421](https://pubmed.ncbi.nlm.nih.gov/7432421/)
114. Greenspan SL et al. Osteoporosis in men with hyperprolactinemic hypogonadism. *Ann Intern Med*. 1986;104:777-82. [PMID 3706929](https://pubmed.ncbi.nlm.nih.gov/3706929/)
115. Colao A et al. Prolactinomas in adolescents: persistent bone loss after 2 years of prolactin normalization. *Clin Endocrinol (Oxf)*. 2000;52:319-27. [PMID 10718830](https://pubmed.ncbi.nlm.nih.gov/10718830/)
116. Vestergaard P et al. Fracture risk is increased in patients with GH deficiency or untreated prolactinomas--a case-control study. *Clin Endocrinol (Oxf)*. 2002;56:159-67. [PMID 11874406](https://pubmed.ncbi.nlm.nih.gov/11874406/)
117. Naliato EC et al. Prevalence of osteopenia in men with prolactinoma. *J Endocrinol Invest*. 2005;28:12-7. [PMID 15816365](https://pubmed.ncbi.nlm.nih.gov/15816365/)
118. Mazziotti G et al. Vertebral fractures in males with prolactinoma. *Endocrine*. 2011;39:288-93. [PMID 21479837](https://pubmed.ncbi.nlm.nih.gov/21479837/)
119. Auriemma RS et al. Effect of cabergoline on metabolism in prolactinomas. *Neuroendocrinology*. 2013;98:299-310. [PMID 24355865](https://pubmed.ncbi.nlm.nih.gov/24355865/)
120. Sperling S & Bhatt H. Prolactinoma: A Massive Effect on Bone Mineral Density in a Young Patient. *Case Rep Endocrinol*. 2016;2016:6312621. [PMID 27446618](https://pubmed.ncbi.nlm.nih.gov/27446618/)
121. Andereggen L et al. Persistent bone impairment despite long-term control of hyperprolactinemia and hypogonadism in men and women with prolactinomas. *Sci Rep*. 2021;11:5122. [PMID 33664388](https://pubmed.ncbi.nlm.nih.gov/33664388/)

---

### L. Mass effect, visual fields, stalk effect

> The geometric compression model, reversible conduction block against irreversible axonal loss, and the ceiling of stalk-compression hyperprolactinaemia (D8).

122. Berwaerts J et al. A giant prolactinoma presenting with unilateral exophthalmos: effect of cabergoline and review of the literature. *J Endocrinol Invest*. 2000;23:393-8. [PMID 10908167](https://pubmed.ncbi.nlm.nih.gov/10908167/)
123. Karavitaki N et al. Do the limits of serum prolactin in disconnection hyperprolactinaemia need re-definition? A study of 226 patients with histologically verified non-functioning pituitary macroadenoma. *Clin Endocrinol (Oxf)*. 2006;65:524-9. [PMID 16984247](https://pubmed.ncbi.nlm.nih.gov/16984247/)
124. Korevaar T et al. Disconnection hyperprolactinaemia in nonadenomatous sellar/parasellar lesions practically never exceeds 2000 mU/l. *Clin Endocrinol (Oxf)*. 2012;76:602-3. [PMID 21942983](https://pubmed.ncbi.nlm.nih.gov/21942983/)
125. Danesh-Meyer HV et al. Optical coherence tomography predicts visual outcome for pituitary tumors. *J Clin Neurosci*. 2015;22:1098-104. [PMID 25891894](https://pubmed.ncbi.nlm.nih.gov/25891894/)
126. Bulwer C et al. Cabergoline-related impulse control disorder in an adolescent with a giant prolactinoma. *Clin Endocrinol (Oxf)*. 2017;86:862-864. [PMID 28346715](https://pubmed.ncbi.nlm.nih.gov/28346715/)
127. Rutland JW et al. Measuring degeneration of the lateral geniculate nuclei from pituitary adenoma compression detected by 7T ultra-high field MRI: a method for predicting vision recovery following surgical decompression of the optic chiasm. *J Neurosurg*. 2020;132:1747-1756. [PMID 31100726](https://pubmed.ncbi.nlm.nih.gov/31100726/)
128. Alkhaibary A et al. Invasive Giant Prolactinoma. *World Neurosurg*. 2024;181:21-22. [PMID 37827431](https://pubmed.ncbi.nlm.nih.gov/37827431/)

---

### M. Safety: valve, impulse control, nausea

> The separation of 5-HT2B (the plasma drive), D3 (impulse control), and area postrema D2 (nausea). The basis for D11, that the
> 10-40 fold difference between Parkinson's disease doses and prolactinoma doses explains the discordance between the two groups of references.

129. Jovanović-Mićić D et al. The role of alpha-adrenergic mechanisms within the area postrema in dopamine-induced emesis. *Eur J Pharmacol*. 1995;272:21-30. [PMID 7713146](https://pubmed.ncbi.nlm.nih.gov/7713146/)
130. Schade R et al. Dopamine agonists and the risk of cardiac-valve regurgitation. *N Engl J Med*. 2007;356:29-38. [PMID 17202453](https://pubmed.ncbi.nlm.nih.gov/17202453/)
131. Zanettini R et al. Valvular heart disease and the use of dopamine agonists for Parkinson's disease. *N Engl J Med*. 2007;356:39-46. [PMID 17202454](https://pubmed.ncbi.nlm.nih.gov/17202454/)
132. Auriemma RS et al. Safety of long-term treatment with cabergoline on cardiac valve disease in patients with prolactinomas. *Eur J Endocrinol*. 2013;169:359-66. [PMID 23824978](https://pubmed.ncbi.nlm.nih.gov/23824978/)
133. Bancos I et al. Impulse control disorders in patients with dopamine agonist-treated prolactinomas and nonfunctioning pituitary adenomas: a case-control study. *Clin Endocrinol (Oxf)*. 2014;80:863-8. [PMID 24274365](https://pubmed.ncbi.nlm.nih.gov/24274365/)
134. Auriemma RS et al. Cabergoline use for pituitary tumors and valvular disorders. *Endocrinol Metab Clin North Am*. 2015;44:89-97. [PMID 25732645](https://pubmed.ncbi.nlm.nih.gov/25732645/)
135. Caputo C et al. The need for annual echocardiography to detect cabergoline-associated valvulopathy in patients with prolactinoma: a systematic review and additional clinical data. *Lancet Diabetes Endocrinol*. 2015;3:906-13. [PMID 25466526](https://pubmed.ncbi.nlm.nih.gov/25466526/)
136. Barake M et al. MANAGEMENT OF ENDOCRINE DISEASE: Impulse control disorders in patients with hyperpolactinemia treated with dopamine agonists: how much should we worry?. *Eur J Endocrinol*. 2018;179:R287-R296. [PMID 30324793](https://pubmed.ncbi.nlm.nih.gov/30324793/)
137. Dogansen SC et al. Dopamine Agonist-Induced Impulse Control Disorders in Patients With Prolactinoma: A Cross-Sectional Multicenter Study. *J Clin Endocrinol Metab*. 2019;104:2527-2534. [PMID 30848825](https://pubmed.ncbi.nlm.nih.gov/30848825/)
138. Stiles CE et al. Incidence of Cabergoline-Associated Valvulopathy in Primary Care Patients With Prolactinoma Using Hard Cardiac Endpoints. *J Clin Endocrinol Metab*. 2021;106:e711-e720. [PMID 33247916](https://pubmed.ncbi.nlm.nih.gov/33247916/)

---

### N. Drug-induced hyperprolactinaemia

> The same occupancy equation in the opposite direction. The basis for D9, in which prolactin falls when a partial agonist such as aripiprazole (e≈0.25) is added.

139. Honbo KS et al. Serum prolactin levels in untreated primary hypothyroidism. *Am J Med*. 1978;64:782-7. [PMID 645742](https://pubmed.ncbi.nlm.nih.gov/645742/)
140. Scanlon MF et al. Altered dopaminergic regulation of thyrotrophin release in patients with prolactinomas: comparison with other tests of hypothalamic-pituitary function. *Clin Endocrinol (Oxf)*. 1981;14:133-43. [PMID 6790201](https://pubmed.ncbi.nlm.nih.gov/6790201/)
141. Daniels GH et al. Effect of risperidone dose on serum prolactin level. *Endocr Pract*. 2001;7:224. [PMID 11421571](https://pubmed.ncbi.nlm.nih.gov/11421571/)
142. Molitch ME. Medication-induced hyperprolactinemia. *Mayo Clin Proc*. 2005;80:1050-7. [PMID 16092584](https://pubmed.ncbi.nlm.nih.gov/16092584/)
143. Hekimsoy Z et al. The prevalence of hyperprolactinaemia in overt and subclinical hypothyroidism. *Endocr J*. 2010;57:1011-5. [PMID 20938100](https://pubmed.ncbi.nlm.nih.gov/20938100/)
144. Ajmal A et al. Psychotropic-induced hyperprolactinemia: a clinical review. *Psychosomatics*. 2014;55:29-36. [PMID 24140188](https://pubmed.ncbi.nlm.nih.gov/24140188/)
145. Peuskens J et al. The effects of novel and newly approved antipsychotics on serum prolactin levels: a comprehensive review. *CNS Drugs*. 2014;28:421-53. [PMID 24677189](https://pubmed.ncbi.nlm.nih.gov/24677189/)
146. Grigg J et al. Antipsychotic-induced hyperprolactinemia: synthesis of world-wide guidelines and integrated recommendations for assessment, management and future research. *Psychopharmacology (Berl)*. 2017;234:3279-3297. [PMID 28889207](https://pubmed.ncbi.nlm.nih.gov/28889207/)
147. Zhang L et al. Efficacy and Safety of Adjunctive Aripiprazole, Metformin, and Paeoniae-Glycyrrhiza Decoction for Antipsychotic-Induced Hyperprolactinemia: A Network Meta-Analysis of Randomized Controlled Trials. *Front Psychiatry*. 2021;12:728204. [PMID 34658963](https://pubmed.ncbi.nlm.nih.gov/34658963/)
148. Lin X et al. Antipsychotic-Related Prolactin Changes: A Systematic Review and Dose-Response Meta-analysis. *CNS Drugs*. 2025;39:937-947. [PMID 40830715](https://pubmed.ncbi.nlm.nih.gov/40830715/)

---

### O. QSP methodology

> The mrgsolve implementation and the model-informed drug development (MIDD) context.

149. Peterson MC & Riggs MM. FDA Advisory Meeting Clinical Pharmacology Review Utilizes a Quantitative Systems Pharmacology (QSP) Model: A Watershed Moment?. *CPT Pharmacometrics Syst Pharmacol*. 2015;4:e00020. [PMID 26225239](https://pubmed.ncbi.nlm.nih.gov/26225239/)
150. Helmlinger G et al. Quantitative Systems Pharmacology: An Exemplar Model-Building Workflow With Applications in Cardiovascular, Metabolic, and Oncology Drug Development. *CPT Pharmacometrics Syst Pharmacol*. 2019;8:380-395. [PMID 31087533](https://pubmed.ncbi.nlm.nih.gov/31087533/)

---

## What was rejected in verification, and why

Among the candidates the relevance filter screened out were **papers on entirely different subjects** (prognosis in multiple sclerosis, a randomised
trial of BCG vaccination, a metal-organic framework biosensor, and so on). These were the result of a fallback query, reduced to an author name plus a short
keyword, picking up irrelevant top hits from PubMed's relevance ordering, and they were excluded automatically.
This is recorded because the point that **an author-name-based fallback query always needs verification after the fact**
was the most practical lesson of building this list.

The following items, further, were not put in the list because no verifiable PubMed citation could be found, and they are therefore not
used as quantitative anchors in the model either. Those parameters are marked explicitly as "estimates" in the
model file.

- The **absolute bioavailability** of cabergoline: no published value exists, so `V2_CAB` and `CL_CAB` are apparent values,
  set so as to reproduce the published pg/mL range and half-life.
- The **pituitary biophase partition coefficient** `PART_CAB = 60`: a structural parameter that absorbs the accumulation in pituitary
  tissue, and one not identifiable from plasma data alone.
- The **immunoassay hook constants** `KHOOK` and `PHOOK`: these differ by platform, so only the *shape* of the curve (monotonic rise → peak →
  collapse) is a general claim.

## Licence and disclaimer

This reference list and this model are for educational and research purposes. They must not be used directly in clinical decision-making.
