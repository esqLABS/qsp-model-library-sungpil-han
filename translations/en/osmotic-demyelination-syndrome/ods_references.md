# Osmotic Demyelination Syndrome — References

> **Model directory:** `osmotic-demyelination-syndrome/` · **Abbreviation:** ODS
> **PMID verification:** every PMID below was looked up with NCBI E-utilities (`esearch`/`esummary`) and
> its author, year, journal, and title checked against the record. No citation here was written from memory.

The `[model]` mark after each entry indicates **which part** of the model that reference supports or refutes.

---

## 1. The original description and the classical pathology

1. Adams RD, Victor M, Mancall EL. Central pontine myelinolysis: a hitherto undescribed disease occurring in alcoholic and malnourished patients. *AMA Arch Neurol Psychiatry* 1959. [PMID 13616772](https://pubmed.ncbi.nlm.nih.gov/13616772/) — the original description. The two risk factors, alcoholism and malnutrition, are already there in the first sentence. `[model: the basis for the existence of FOSM and FNUT]`
2. Wright DG, Laureno R, Victor M. Pontine and extrapontine myelinolysis. *Brain* 1979. [PMID 455045](https://pubmed.ncbi.nlm.nih.gov/455045/) — a systematic description of the extrapontine lesions. `[model: WEXP = 0.55]`
3. Norenberg MD, Leslie KO, Robertson AS. Association between rise in serum sodium and central pontine myelinolysis. *Ann Neurol* 1982. [PMID 7073246](https://pubmed.ncbi.nlm.nih.gov/7073246/) — the first clinical link to "rate" as the cause.
4. Kleinschmidt-DeMasters BK, Norenberg MD. Rapid correction of hyponatremia causes demyelination: relation to central pontine myelinolysis. *Science* 1981. [PMID 7466381](https://pubmed.ncbi.nlm.nih.gov/7466381/) — causation established in the rat model. The experiment this whole model sets out to reproduce.
5. Sterns RH, Riggs JE, Schochet SS. Osmotic demyelination syndrome following correction of hyponatremia. *N Engl J Med* 1986. [PMID 3713747](https://pubmed.ncbi.nlm.nih.gov/3713747/) — a clinical case series. The origin of the name "osmotic demyelination syndrome".
6. Karp BI, Laureno R. Pontine and extrapontine myelinolysis: a neurologic disorder following rapid correction of hyponatremia. *Medicine (Baltimore)* 1993. [PMID 8231786](https://pubmed.ncbi.nlm.nih.gov/8231786/)
7. Martin RJ. Central pontine and extrapontine myelinolysis: the osmotic demyelination syndromes. *J Neurol Neurosurg Psychiatry* 2004. [PMID 15316041](https://pubmed.ncbi.nlm.nih.gov/15316041/)
8. Brown WD. Osmotic demyelination disorders: central pontine and extrapontine myelinolysis. *Curr Opin Neurol* 2000. [PMID 11148672](https://pubmed.ncbi.nlm.nih.gov/11148672/)
9. King JD, Rosner MH. Osmotic demyelination syndrome. *Am J Med Sci* 2010. [PMID 20453633](https://pubmed.ncbi.nlm.nih.gov/20453633/)
10. Murase T, Sugimura Y, Takefuji S, Oiso Y, Murata Y. Mechanisms and therapy of osmotic demyelination. *Am J Med* 2006. [PMID 16843088](https://pubmed.ncbi.nlm.nih.gov/16843088/)

---

## 2. Brain osmolyte adaptation — the core of this model

11. Verbalis JG, Gullans SR. Hyponatremia causes large sustained reductions in brain content of multiple organic osmolytes in rats. *Brain Res* 1991. [PMID 1817731](https://pubmed.ncbi.nlm.nih.gov/1817731/) — **the basis for calibrating β_i in the model.** In chronic hyponatraemia the organic osmolytes are lost in bulk and the loss is sustained. At [Na] 110 the model gives total organic osmolytes −40% and myo-inositol −66%.
12. Verbalis JG, Gullans SR. Rapid correction of hyponatremia produces differential effects on brain osmolyte and electrolyte reaccumulation in rats. *Brain Res* 1993. [PMID 8096428](https://pubmed.ncbi.nlm.nih.gov/8096428/) — **the paper that is the reason this model exists.** The electrolytes recover within 24 hours, the organic osmolytes over about 5 days. `[model: TAUEFF 0.35 d vs TAUINF 1.25 d × the transcriptional delay]`
13. Lien YH, Shapiro JI, Chan L. Study of brain electrolytes and organic osmolytes during correction of chronic hyponatremia. Implications for the pathogenesis of central pontine myelinolysis. *J Clin Invest* 1991. [PMID 2056123](https://pubmed.ncbi.nlm.nih.gov/2056123/) — the decisive study, measuring the inorganic ions and the organic osmolytes separately and proposing the mechanism of CPM.
14. Verbalis JG, Drutarosky MD. Adaptation to chronic hypoosmolality in rats. *Kidney Int* 1988. [PMID 3172643](https://pubmed.ncbi.nlm.nih.gov/3172643/) — the observation that the water content of the adapted brain is almost normal. `[model: after adaptation at [Na] 110, brain water 80.08 against a normal 80.00]`
15. Lien YH, Shapiro JI, Chan L. Effects of hypernatremia on organic brain osmoles. *J Clin Invest* 1990. [PMID 2332498](https://pubmed.ncbi.nlm.nih.gov/2332498/) — osmolyte kinetics in adaptation in the opposite direction (hypernatraemia).
16. Videen JS, Michaelis T, Pinto P, Ross BD. Human cerebral osmolytes during chronic hyponatremia. A proton magnetic resonance spectroscopy study. *J Clin Invest* 1995. [PMID 7860762](https://pubmed.ncbi.nlm.nih.gov/7860762/) — **the only series measured directly in humans.** The model's Ω has never been measured in a human, and MRS myo-inositol is one of its terms.
17. Häussinger D, Laubenberger J, vom Dahl S, et al. Proton magnetic resonance spectroscopy studies on human brain myo-inositol in hypo-osmolarity and hepatic encephalopathy. *Gastroenterology* 1994. [PMID 7926510](https://pubmed.ncbi.nlm.nih.gov/7926510/) — **the basis for the liver disease phenotype.** In hepatic encephalopathy myo-inositol is already depleted. `[model: cirrhosis phenotype INS0 7.0 → 4.2]`
18. Restuccia T, Gómez-Ansón B, Guevara M, et al. Effects of dilutional hyponatremia on brain organic osmolytes and water content in patients with cirrhosis. *Hepatology* 2004. [PMID 15185302](https://pubmed.ncbi.nlm.nih.gov/15185302/)
19. Massieu L, Montiel T, Robles G, Quesada O. Brain amino acids during hyponatremia in vivo: clinical observations and experimental studies. *Neurochem Res* 2004. [PMID 14992265](https://pubmed.ncbi.nlm.nih.gov/14992265/)
20. Sterns RH, Silver SM. Brain volume regulation in response to hypo-osmolality and its correction. *Am J Med* 2006. [PMID 16843080](https://pubmed.ncbi.nlm.nih.gov/16843080/)
21. Gankam Kengne F, Decaux G. Hyponatremia and the Brain. *Kidney Int Rep* 2018. [PMID 29340311](https://pubmed.ncbi.nlm.nih.gov/29340311/)
22. Gankam Kengne F. Adaptation of the Brain to Hyponatremia and Its Clinical Implications. *J Clin Med* 2023. [PMID 36902500](https://pubmed.ncbi.nlm.nih.gov/36902500/)

---

## 3. The molecular asymmetry of osmolyte transport — channel against carrier

> This section carries the structural argument of the model: **what takes them out is a channel; what brings them back in is a carrier that has to be transcribed.**

23. Qiu Z, Dubin AE, Mathur J, et al. SWELL1, a plasma membrane protein, is an essential component of volume-regulated anion channel. *Cell* 2014. [PMID 24725410](https://pubmed.ncbi.nlm.nih.gov/24725410/) — the molecular identity of VRAC (LRRC8A). `[model: the efflux arm, TAUEFF]`
24. Voss FK, Ullrich F, Münch J, et al. Identification of LRRC8 heteromers as an essential component of the volume-regulated anion channel VRAC. *Science* 2014. [PMID 24790029](https://pubmed.ncbi.nlm.nih.gov/24790029/)
25. Thöne FMB, et al. LRRC8/VRAC chloride and metabolite channels in signaling and volume regulation. *Trends Biochem Sci* 2025. [PMID 40769835](https://pubmed.ncbi.nlm.nih.gov/40769835/) — the most recent review.
26. Miyakawa H, Woo SK, Dahl SC, Handler JS, Kwon HM. Tonicity-responsive enhancer binding protein, a rel-like protein that stimulates transcription in response to hypertonicity. *Proc Natl Acad Sci U S A* 1999. [PMID 10051678](https://pubmed.ncbi.nlm.nih.gov/10051678/) — **TonEBP/NFAT5.** The basis for locating the reason reaccumulation is slow in transcription. `[model: TAUTON 0.25 d → TAUSMIT 1.0 d]`
27. Ko BC, Lam AK, Kapus A, Fan L, Chung SK, Chung SS. Fyn and p38 signaling are both required for maximal hypertonic activation of the osmotic response element-binding protein/tonicity-responsive enhancer-binding protein. *J Biol Chem* 2002. [PMID 12359721](https://pubmed.ncbi.nlm.nih.gov/12359721/)
28. Burg MB, Ferraris JD, Dmitrieva NI. Cellular response to hyperosmotic stresses. *Physiol Rev* 2007. [PMID 17928589](https://pubmed.ncbi.nlm.nih.gov/17928589/) — the standard review of the TonEBP → SMIT1/TauT/BGT1/AR target gene set.
29. Berry GT, Mallee JJ, Kwon HM, et al. The human osmoregulatory Na+/myo-inositol cotransporter gene (SLC5A3): molecular cloning and localization to chromosome 21. *Genomics* 1995. [PMID 7789985](https://pubmed.ncbi.nlm.nih.gov/7789985/) — the cloning of SMIT1. That it is **Na⁺-coupled** is the molecular basis for the risk carried by hypokalaemia.
30. Bissonnette P, Lahjouji K, Coady MJ, Lapointe JY. Effects of hyperosmolarity on the Na+-myo-inositol cotransporter SMIT2 stably transfected in the Madin-Darby canine kidney cell line. *Am J Physiol Cell Physiol* 2008. [PMID 18650262](https://pubmed.ncbi.nlm.nih.gov/18650262/)
31. Pasantes-Morales H. Taurine Homeostasis and Volume Control. *Adv Neurobiol* 2017. [PMID 28828605](https://pubmed.ncbi.nlm.nih.gov/28828605/) — the role of TauT and taurine in volume regulation.
32. Morán J, Maar T, Pasantes-Morales H. Impaired cell volume regulation in taurine deficient cultured astrocytes. *Neurochem Res* 1994. [PMID 8065498](https://pubmed.ncbi.nlm.nih.gov/8065498/) — direct evidence that volume regulation fails when the osmolyte pool is empty.
33. Bourque CW. Central mechanisms of osmosensation and systemic osmoregulation. *Nat Rev Neurosci* 2008. [PMID 18509340](https://pubmed.ncbi.nlm.nih.gov/18509340/) — the osmoreceptors (OVLT/SFO). `[model: OSMTHR 280, GAINOSM 0.45]`
34. Papadopoulos MC, Verkman AS. Aquaporin water channels in the nervous system. *Nat Rev Neurosci* 2013. [PMID 23481483](https://pubmed.ncbi.nlm.nih.gov/23481483/) — AQP4. The basis for the assumption that water equilibrates within seconds (brain water is computed algebraically in the model).

---

## 4. The order of cellular injury — the astrocyte first

35. Gankam Kengne F, Nicaise C, Soupart A, et al. Astrocytes are an early target in osmotic demyelination syndrome. *J Am Soc Nephrol* — for the related series of studies see entries 36–39. (The basis for the `AST → OLI` order in this model)
36. Gankam-Kengne F, Soupart A, Pochet R, Brion JP, Decaux G. Minocycline protects against neurologic complications of rapid correction of hyponatremia. *J Am Soc Nephrol* 2010. [PMID 21051736](https://pubmed.ncbi.nlm.nih.gov/21051736/) — microglial inhibition is protective. `[model: MINOON, EMAXMIN]`
37. Takefuji S, Murase T, Sugimura Y, et al. Role of microglia in the pathogenesis of osmotic-induced demyelination. *Exp Neurol* 2007. [PMID 17126835](https://pubmed.ncbi.nlm.nih.gov/17126835/)
38. Iwama S, Sugimura Y, Suzuki H, et al. Time-dependent changes in proinflammatory and neurotrophic responses of microglia and astrocytes in a rat model of osmotic demyelination syndrome. *Glia* 2011. [PMID 21264951](https://pubmed.ncbi.nlm.nih.gov/21264951/) — **the direct basis for the model's time order (astrocyte → microglia → demyelination).**
39. Sugimura Y, Murase T, Takefuji S, et al. Protective effect of dexamethasone on osmotic-induced demyelination in rats. *Exp Neurol* 2005. [PMID 15698632](https://pubmed.ncbi.nlm.nih.gov/15698632/) — `[model: DEXON, EMAXDEX — acting on the BBB arm only]`

---

## 5. Urea and relowering — the animal evidence and the model's honest negative result

40. Soupart A, Penninckx R, Stenuit A, Decaux G. Azotemia (48 h) decreases the risk of brain damage in rats after correction of chronic hyponatremia. *Brain Res* 2000. [PMID 10661508](https://pubmed.ncbi.nlm.nih.gov/10661508/)
41. Gankam Kengne F, Couturier BS, Soupart A, Decaux G. Urea minimizes brain complications following rapid correction of chronic hyponatremia compared with vasopressin antagonist or hypertonic saline. *Kidney Int* 2015. [PMID 25100046](https://pubmed.ncbi.nlm.nih.gov/25100046/) — **the most exposed assumption in the model.** On a pure effective-osmolality accounting, urea crosses the barrier and so cancels on both sides and cannot protect. The model therefore lets urea act **on the BBB and microglial arms only** (where those studies actually measured). If urea is shown to reduce the osmotic injury itself, this model is wrong.
42. Soupart A, Ngassa M, Decaux G. Reinduction of hyponatremia improves survival in rats with myelinolysis-related neurologic symptoms. *J Neuropathol Exp Neurol* 1996. [PMID 8627349](https://pubmed.ncbi.nlm.nih.gov/8627349/) — the origin of relowering.
43. Gankam Kengne F, Soupart A, Pochet R, Brion JP, Decaux G. Re-induction of hyponatremia after rapid overcorrection of hyponatremia reduces mortality in rats. *Kidney Int* 2009. [PMID 19606078](https://pubmed.ncbi.nlm.nih.gov/19606078/) — `[model: the §6 relowering deadline experiment]`
44. Soupart A, Penninckx R, Namias B, et al. — for the related clinical relowering cases, read alongside the desmopressin series at entries 45–47.
45. Decaux G, Andres C, Gankam Kengne F, Soupart A. Treatment of euvolemic hyponatremia in the intensive care unit by urea. *Crit Care* 2010. [PMID 20946646](https://pubmed.ncbi.nlm.nih.gov/20946646/)
46. Decaux G, Genette F. Urea for long-term treatment of syndrome of inappropriate secretion of antidiuretic hormone. *Br Med J (Clin Res Ed)* 1981. [PMID 6794768](https://pubmed.ncbi.nlm.nih.gov/6794768/)
47. Lockett J, Berkman KE, Dimeski G, Russell AW, Inder WJ. Urea treatment in fluid restriction-refractory hyponatraemia. *Clin Endocrinol (Oxf)* 2019. [PMID 30614552](https://pubmed.ncbi.nlm.nih.gov/30614552/)

---

## 6. Desmopressin — taking the kidney out of the loop

> What §4 and §5 of the model claim: **what sets the correction rate is not the prescription but the kidney.**
> Treatment is therefore not a matter of raising the sodium but of stopping the kidney.

48. Perianayagam A, Sterns RH, Silver SM, et al. DDAVP is effective in preventing and reversing inadvertent overcorrection of hyponatremia. *Clin J Am Soc Nephrol* 2008. [PMID 18235152](https://pubmed.ncbi.nlm.nih.gov/18235152/)
49. Sood L, Sterns RH, Hix JK, Silver SM, Chen L. Hypertonic saline and desmopressin: a simple strategy for safe correction of severe hyponatremia. *Am J Kidney Dis* 2013. [PMID 23266328](https://pubmed.ncbi.nlm.nih.gov/23266328/) — **the origin of the proactive DDAVP clamp.** `[model: S08]`
50. Rafat C, Schortgen F, Gaudry S, et al. Use of desmopressin acetate in severe hyponatremia in the intensive care unit. *Clin J Am Soc Nephrol* 2014. [PMID 24262506](https://pubmed.ncbi.nlm.nih.gov/24262506/)
51. MacMillan TE, Cavalcanti RB. Outcomes in Severe Hyponatremia Treated With and Without Desmopressin. *Am J Med* 2018. [PMID 29061503](https://pubmed.ncbi.nlm.nih.gov/29061503/)
52. Sterns RH. Treatment of Severe Hyponatremia. *Clin J Am Soc Nephrol* 2018. [PMID 29295830](https://pubmed.ncbi.nlm.nih.gov/29295830/)
53. Sterns RH. Adverse Consequences of Overly-Rapid Correction of Hyponatremia. *Front Horm Res* 2019. [PMID 32097948](https://pubmed.ncbi.nlm.nih.gov/32097948/)

---

## 7. Correction rate, overcorrection, and their consequences

54. Sterns RH, Cappuccio JD, Silver SM, Cohen EP. Neurologic sequelae after treatment of severe hyponatremia: a multicenter perspective. *J Am Soc Nephrol* 1994. [PMID 8025225](https://pubmed.ncbi.nlm.nih.gov/8025225/) — the clinical basis for the rate limit.
55. Ayus JC, Krothapalli RK, Arieff AI. Treatment of symptomatic hyponatremia and its relation to brain damage. A prospective study. *N Engl J Med* 1987. [PMID 3309659](https://pubmed.ncbi.nlm.nih.gov/3309659/)
56. Ayus JC, Krothapalli RK, Armstrong DL. Changing concepts in treatment of severe symptomatic hyponatremia. Rapid correction and possible relation to central pontine myelinolysis. *Am J Med* 1985. [PMID 4014266](https://pubmed.ncbi.nlm.nih.gov/4014266/)
57. Berl T. Treating hyponatremia: damned if we do and damned if we don't. *Kidney Int* 1990. [PMID 2179612](https://pubmed.ncbi.nlm.nih.gov/2179612/) — the classic formulation of the dilemma. This model shows that the dilemma points in opposite directions in acute and in chronic disease (§12).
58. George JC, Zafar W, Bucaloiu ID, Chang AR. Risk Factors and Outcomes of Rapid Correction of Severe Hyponatremia. *Clin J Am Soc Nephrol* 2018. [PMID 29871886](https://pubmed.ncbi.nlm.nih.gov/29871886/)
59. Woodfine JD, Sood MM, MacMillan TE, Cavalcanti RB, van Walraven C. Criteria for Hyponatremic Overcorrection: Systematic Review and Cohort Study of Emergently Ill Patients. *J Gen Intern Med* 2020. [PMID 31452039](https://pubmed.ncbi.nlm.nih.gov/31452039/)
60. Mohmand HK, Issa D, Ahmad Z, Cappuccio JD, Kouides RW, Sterns RH. Hypertonic saline for hyponatremia: risk of inadvertent overcorrection. *Clin J Am Soc Nephrol* 2007. [PMID 17913972](https://pubmed.ncbi.nlm.nih.gov/17913972/)
61. Baek SH, Jo YH, Ahn S, et al. Risk of Overcorrection in Rapid Intermittent Bolus vs Slow Continuous Infusion Therapies of Hypertonic Saline for Patients With Symptomatic Hyponatremia: The SALSA Randomized Clinical Trial. *JAMA Intern Med* 2021. [PMID 33104189](https://pubmed.ncbi.nlm.nih.gov/33104189/)
62. Yun G, et al. Risk factors for undercorrection in severe hyponatremia: a post-hoc analysis of the SALSA trial. *Kidney Res Clin Pract* 2026. [PMID 42311099](https://pubmed.ncbi.nlm.nih.gov/42311099/)
63. Refardt J, et al. A Randomized Trial of Targeted Hyponatremia Correction in Hospitalized Patients. *NEJM Evid* 2026. [PMID 41733398](https://pubmed.ncbi.nlm.nih.gov/41733398/)
64. Tandukar S, Sterns RH, Rondon-Berrios H. Osmotic Demyelination Syndrome following Correction of Hyponatremia by ≤10 mEq/L per Day. *Kidney360* 2021. [PMID 35373113](https://pubmed.ncbi.nlm.nih.gov/35373113/) — **exactly the group of patients §11 of the model predicts.** Cases in which demyelination occurred even though the prescription was adhered to.
65. Sterns RH. Treatment of hyponatremia. *Curr Opin Nephrol Hypertens* 2010. [PMID 20539224](https://pubmed.ncbi.nlm.nih.gov/20539224/)
66. Achinger SG, Ayus JC. Treatment of Hyponatremic Encephalopathy in the Critically Ill. *Crit Care Med* 2017. [PMID 28704229](https://pubmed.ncbi.nlm.nih.gov/28704229/)
67. Moritz ML, Ayus JC. 100 cc 3% sodium chloride bolus: a novel treatment for hyponatremic encephalopathy. *Metab Brain Dis* 2010. [PMID 20221678](https://pubmed.ncbi.nlm.nih.gov/20221678/) — `[model: the §14 Adrogué–Madias cross-check]`

---

## 8. Guidelines and reviews

68. Spasovski G, Vanholder R, Allolio B, et al. Clinical practice guideline on diagnosis and treatment of hyponatraemia. *Nephrol Dial Transplant* 2014. [PMID 24569496](https://pubmed.ncbi.nlm.nih.gov/24569496/) — the European guideline (a 150 mL bolus of 3% NaCl × 3, targeting a rise of 5 mmol/L). `[model: the object of the §14 cross-check]`
69. Verbalis JG, Goldsmith SR, Greenberg A, et al. Diagnosis, evaluation, and treatment of hyponatremia: expert panel recommendations. *Am J Med* 2013. [PMID 24074529](https://pubmed.ncbi.nlm.nih.gov/24074529/) — the American expert panel. Normal risk ≤10–12, high risk ≤8 mmol/L/24 h. **The model derives these two numbers from a single threshold (Ω* = 8) and a single carrier (FOSM).**
70. Hoorn EJ, Zietse R. Diagnosis and Treatment of Hyponatremia: Compilation of the Guidelines. *J Am Soc Nephrol* 2017. [PMID 28174217](https://pubmed.ncbi.nlm.nih.gov/28174217/)
71. Sterns RH. Disorders of plasma sodium — causes, consequences, and correction. *N Engl J Med* 2015. [PMID 25551526](https://pubmed.ncbi.nlm.nih.gov/25551526/)
72. Adrogué HJ, Madias NE. Hyponatremia. *N Engl J Med* 2000. [PMID 10824078](https://pubmed.ncbi.nlm.nih.gov/10824078/) — the origin of the Adrogué–Madias formula.
73. Seay NW, Lehrich RW, Greenberg A. Diagnosis and Management of Disorders of Body Tonicity — Hyponatremia and Hypernatremia: Core Curriculum 2020. *Am J Kidney Dis* 2020. [PMID 31606238](https://pubmed.ncbi.nlm.nih.gov/31606238/)
74. Sterns RH, Nigwekar SU, Hix JK. The treatment of hyponatremia. *Semin Nephrol* 2009. [PMID 19523575](https://pubmed.ncbi.nlm.nih.gov/19523575/)
75. Overgaard-Steensen C, Ring T. Clinical review: practical approach to hyponatraemia and hypernatraemia in critically ill patients. *Crit Care* 2013. [PMID 23672688](https://pubmed.ncbi.nlm.nih.gov/23672688/)
76. Ellison DH, Berl T. Clinical practice. The syndrome of inappropriate antidiuresis. *N Engl J Med* 2007. [PMID 17507705](https://pubmed.ncbi.nlm.nih.gov/17507705/)
77. Rosner MH, et al. Syndrome of Inappropriate Antidiuresis. *J Am Soc Nephrol* 2025. [PMID 39621420](https://pubmed.ncbi.nlm.nih.gov/39621420/)
78. Prince R, et al. The treatment of acute symptomatic hyponatraemia in the hospital setting. *Best Pract Res Clin Endocrinol Metab* 2026. [PMID 41402221](https://pubmed.ncbi.nlm.nih.gov/41402221/)

---

## 9. The quantitative water-and-salt layer

79. Edelman IS, Leibman J, O'Meara MP, Birkenfeld LW. Interrelations between serum sodium concentration, serum osmolarity and total exchangeable sodium, total exchangeable potassium and total body water. *J Clin Invest* 1958. [PMID 13575523](https://pubmed.ncbi.nlm.nih.gov/13575523/) — **the model's [Na]s = 1.11(Na_e+K_e)/TBW − 25.6.** That potassium is in the numerator is the whole of §7.
80. Nguyen MK, Kurtz I. New insights into the pathophysiology of the dysnatremias: a quantitative analysis. *Am J Physiol Renal Physiol* 2004. [PMID 15271684](https://pubmed.ncbi.nlm.nih.gov/15271684/)
81. Kurtz I, Nguyen MK. Evolving concepts in the quantitative analysis of the determinants of the plasma water sodium concentration. *Kidney Int* 2005. [PMID 16221198](https://pubmed.ncbi.nlm.nih.gov/16221198/) — the osmotically inactive sodium store. `[model: FECF 0.65]`
82. Rose BD. New approach to disturbances in the plasma sodium concentration. *Am J Med* 1986. [PMID 3799631](https://pubmed.ncbi.nlm.nih.gov/3799631/)
83. Barsoum NR, Levine BS. Current prescriptions for the correction of hyponatraemia and hypernatraemia: are they too simple? *Nephrol Dial Transplant* 2002. [PMID 12105238](https://pubmed.ncbi.nlm.nih.gov/12105238/) — the point that the formulae ignore the urine. §4 of the model shows that this is precisely the dominant term.
84. Sterns RH. Formulas for fixing serum sodium: curb your enthusiasm. *Clin Kidney J* 2016. [PMID 27478590](https://pubmed.ncbi.nlm.nih.gov/27478590/)
85. Furst H, Hallows KR, Post J, Chen S, Kotzker W, Goldfarb S, Ziyadeh F, Neilson EG. The urine/plasma electrolyte ratio: a predictive guide to water restriction. *Am J Med Sci* 2000. [PMID 10768609](https://pubmed.ncbi.nlm.nih.gov/10768609/) — **the rule the model derives.** If (U_Na+U_K)/S_Na > 1, water restriction will fail. The model's SIADH patient has a ratio of 2.0 at presentation.
86. Berl T. Impact of solute intake on urine flow and water excretion. *J Am Soc Nephrol* 2008. [PMID 18337482](https://pubmed.ncbi.nlm.nih.gov/18337482/) — solute intake governs free water excretion. `[model: beer potomania phenotype PUREA 110 → 400]`
87. Robertson GL. Thirst and vasopressin function in normal and disordered states of water balance. *J Lab Clin Med* 1983. [PMID 6338137](https://pubmed.ncbi.nlm.nih.gov/6338137/) — `[model: OSMTHIRST 292, KTHIRST]`

---

## 10. Phenotypes and risk factors

88. Lohr JW. Osmotic demyelination syndrome following correction of hyponatremia: association with hypokalemia. *Am J Med* 1994. [PMID 8192171](https://pubmed.ncbi.nlm.nih.gov/8192171/) — **in the model hypokalaemia is dangerous by two routes:** it raises the sodium through the Edelman numerator (§7), and it lowers the driving force of the Na⁺-coupled osmolyte carriers.
89. Lee EM, Kang JK, Yun SC, et al. Risk factors for central pontine and extrapontine myelinolysis following orthotopic liver transplantation. *Eur Neurol* 2009. [PMID 19797900](https://pubmed.ncbi.nlm.nih.gov/19797900/)
90. Yun BC, Kim WR, Benson JT, et al. Impact of pretransplant hyponatremia on outcome following liver transplantation. *Hepatology* 2009. [PMID 19402063](https://pubmed.ncbi.nlm.nih.gov/19402063/)
91. Rondon-Berrios H, Velez JCQ. Hyponatremia in Cirrhosis. *Clin Liver Dis* 2022. [PMID 35487602](https://pubmed.ncbi.nlm.nih.gov/35487602/)
92. Thaler SM, Teitelbaum I, Berl T. "Beer potomania" in non-beer drinkers: effect of low dietary solute intake. *Am J Kidney Dis* 1998. [PMID 9631849](https://pubmed.ncbi.nlm.nih.gov/9631849/)
93. Sanghvi SR, Kellerman PS, Nanovic L. Beer potomania: an unusual cause of hyponatremia at high risk of complications from rapid correction. *Am J Kidney Dis* 2007. [PMID 17900468](https://pubmed.ncbi.nlm.nih.gov/17900468/) — `[model: S17 — autonomous overcorrection from the restoration of solute intake alone]`
94. Burst V, Grundmann F, Kubacki T, et al. Thiazide-Associated Hyponatremia, Report of the Hyponatremia Registry. *Am J Nephrol* 2017. [PMID 28419981](https://pubmed.ncbi.nlm.nih.gov/28419981/) — `[model: the THIAZ flag — a raised minimum urine osmolality]`
95. Sailer CO, Winzeler B, Nigro N, et al. Characteristics and outcomes of patients with profound hyponatraemia due to primary polydipsia. *Clin Endocrinol (Oxf)* 2017. [PMID 28556237](https://pubmed.ncbi.nlm.nih.gov/28556237/)
96. Almond CS, Shin AY, Fortescue EB, et al. Hyponatremia among runners in the Boston Marathon. *N Engl J Med* 2005. [PMID 15829535](https://pubmed.ncbi.nlm.nih.gov/15829535/) — the acute phenotype. `[model: §12 — the same sodium, a disease in the opposite direction]`
97. Hew-Butler T, Rosner MH, Fowkes-Godek S, et al. Statement of the Third International Exercise-Associated Hyponatremia Consensus Development Conference. *Clin J Sport Med* 2015. [PMID 26102445](https://pubmed.ncbi.nlm.nih.gov/26102445/)
98. Ayus JC, Arieff AI, Moritz ML. Exercise-associated hyponatremia masquerading as acute mountain sickness. *Clin J Sport Med* 2008. [PMID 18806543](https://pubmed.ncbi.nlm.nih.gov/18806543/)

---

## 11. Course · imaging · outcome

99. Graff-Radford J, Fugate JE, Kaufmann TJ, Mandrekar JN, Rabinstein AA. Clinical and radiologic correlations of central pontine myelinolysis syndrome. *Mayo Clin Proc* 2011. [PMID 21997578](https://pubmed.ncbi.nlm.nih.gov/21997578/)
100. Singh TD, Fugate JE, Rabinstein AA. Central pontine and extrapontine myelinolysis: a systematic review. *Eur J Neurol* 2014. [PMID 25220878](https://pubmed.ncbi.nlm.nih.gov/25220878/)
101. Kallakatta RN, Radhakrishnan A, Fayaz RK, Unnikrishnan JP, Kesavadas C, Sarma SP. Clinical and functional outcome and factors predicting prognosis in osmotic demyelination syndrome. *J Neurol Neurosurg Psychiatry* 2011. [PMID 20826870](https://pubmed.ncbi.nlm.nih.gov/20826870/)
102. Menger H, Jörg J. Outcome of central pontine and extrapontine myelinolysis (n = 44). *J Neurol* 1999. [PMID 10460448](https://pubmed.ncbi.nlm.nih.gov/10460448/) — a substantial proportion recover. `[model: the KMYE and KOPCP remyelination arms]`
103. Louis G, Megarbane B, Lavoué S, et al. Long-term outcome of patients hospitalized in intensive care units with central or extrapontine myelinolysis. *Crit Care Med* 2012. [PMID 22036854](https://pubmed.ncbi.nlm.nih.gov/22036854/)
104. Ruzek KA, Campeau NG, Miller GM. Early diagnosis of central pontine myelinolysis with diffusion-weighted imaging. *AJNR Am J Neuroradiol* 2004. [PMID 14970019](https://pubmed.ncbi.nlm.nih.gov/14970019/) — **the basis for imaging being late.** `[model: §3 — the symptoms come before the MRI]`
105. Chu K, Kang DW, Ko SB, Kim M. Diffusion-weighted MR findings of central pontine and extrapontine myelinolysis. *Acta Neurol Scand* 2001. [PMID 11903095](https://pubmed.ncbi.nlm.nih.gov/11903095/)
106. Laureno R, Illowsky Karp B. Sequential MRI in pontine and extrapontine myelinolysis following rapid correction of hyponatremia. *BMC Res Notes* 2018. [PMID 30290836](https://pubmed.ncbi.nlm.nih.gov/30290836/) — shows the time lag directly on serial MRI.
107. Aegisdottir H, Cooray C, Wirdefeldt K, Piehl F, Sveinsson O. Incidence of osmotic demyelination syndrome in Sweden: A nationwide study. *Acta Neurol Scand* 2019. [PMID 31343728](https://pubmed.ncbi.nlm.nih.gov/31343728/)
108. Manoharan K, et al. Osmotic Demyelination Syndrome Secondary to Hypernatremia and Hypokalemia: A Case Report of Complete Recovery. *J Assoc Physicians India* 2026. [PMID 42543972](https://pubmed.ncbi.nlm.nih.gov/42543972/) — that it occurs even without hyponatraemia (natural enough in the model, where Ω alone is the cause).
109. Yilmaz O, et al. Isolated extrapontine myelinolysis of osmotic demyelination syndrome. *Prague Med Rep* 2013. [PMID 23547724](https://pubmed.ncbi.nlm.nih.gov/23547724/)
110. Vakharia JD, et al. Extrapontine Myelinolysis in a Child with Salt Intoxication. *J Pediatr* 2017. [PMID 28427777](https://pubmed.ncbi.nlm.nih.gov/28427777/)

---

## 12. Vasopressin receptor antagonists (vaptans — the drug that opens the aquaresis)

111. Schrier RW, Gross P, Gheorghiade M, et al. Tolvaptan, a selective oral vasopressin V2-receptor antagonist, for hyponatremia. *N Engl J Med* 2006. [PMID 17105757](https://pubmed.ncbi.nlm.nih.gov/17105757/) — SALT-1/2. `[model: KITLV, S13]`
112. Verbalis JG, Adler S, Schrier RW, et al. Efficacy and safety of oral tolvaptan therapy in patients with the syndrome of inappropriate antidiuretic hormone secretion. *Eur J Endocrinol* 2011. [PMID 21317283](https://pubmed.ncbi.nlm.nih.gov/21317283/)
113. Rondon-Berrios H, Berl T. Vasopressin Receptor Antagonists in Hyponatremia: Uses and Misuses. *Front Med (Lausanne)* 2017. [PMID 28879182](https://pubmed.ncbi.nlm.nih.gov/28879182/)
114. Malhotra I, Gopinath S, Janga KC, et al. Unpredictable nature of tolvaptan in treatment of hypervolemic hyponatremia: case review on role of vaptans. *Case Rep Endocrinol* 2014. [PMID 24511399](https://pubmed.ncbi.nlm.nih.gov/24511399/) — overcorrection provoked by a vaptan.

---

## 13. Epidemiology and prognostic background

115. Corona G, Giuliani C, Parenti G, et al. Moderate hyponatremia is associated with increased risk of mortality: evidence from a meta-analysis. *PLoS One* 2013. [PMID 24367479](https://pubmed.ncbi.nlm.nih.gov/24367479/)
116. Vandergheynst F, Sakr Y, Felleiter P, et al. Incidence and prognosis of dysnatraemia in critically ill patients. *Eur J Clin Invest* 2013. [PMID 23869476](https://pubmed.ncbi.nlm.nih.gov/23869476/)
117. Nzerue CM, Baffoe-Bonnie H, You W, Falana B, Dai S. Predictors of outcome in hospitalized patients with severe hyponatremia. *J Natl Med Assoc* 2003. [PMID 12793790](https://pubmed.ncbi.nlm.nih.gov/12793790/)
118. Ioannou P, et al. Increased Mortality in Elderly Patients Admitted with Hyponatremia. *J Clin Med* 2021. [PMID 34300225](https://pubmed.ncbi.nlm.nih.gov/34300225/)
119. Portales-Castillo I, Sterns RH. Allostasis and the Clinical Manifestations of Mild to Moderate Chronic Hyponatremia: No Good Adaptation Goes Unpunished. *Am J Kidney Dis* 2019. [PMID 30554800](https://pubmed.ncbi.nlm.nih.gov/30554800/) — the view that the adaptation itself carries a price.
120. Ayus JC, Arieff AI. Abnormalities of water metabolism in the elderly. *Semin Nephrol* 1996. [PMID 8829266](https://pubmed.ncbi.nlm.nih.gov/8829266/)

---

## 14. Methodology — QSP and mrgsolve

121. mrgsolve: Simulate from ODE-Based Population PK/PD and QSP Models. <https://mrgsolve.org/>
122. gPKPDviz — an mrgsolve-based Shiny tool for PK/PD simulation. Paper: <https://pmc.ncbi.nlm.nih.gov/articles/PMC10941578/> · code: <https://github.com/Genentech/gPKPDviz/>
123. Vantage Research. QSP in R. <https://vantage-research.net/qsp-in-r/>

---

## Appendix A — The model's falsifiable predictions

Things this model produced **without** taking them from the literature. Each is written in a form that can be wrong.

| # | Prediction | How to refute it |
|---|------|------------------|
| P1 | The two guideline limits (normal risk ≤10–12, high risk ≤8) are not two different rules but the product of one astrocytic threshold and the organic osmolyte carrier capacity (FOSM). | Refuted by showing that the 8 mmol/L/24 h limit is needed in a group of patients whose carrier capacity is normal. |
| P2 | **MRS myo-inositol deficiency at 24–48 hours into correction** separates those who demyelinate better than the 24-hour rise in sodium does. | A prospective cohort with serial ¹H-MRS. The most direct route to testing this model. |
| P3 | In hypovolaemic hyponatraemia, most of the overcorrection that occurs on isotonic saline alone is **autonomous water diuresis from the removal of the AVP stimulus**, not the sodium given. | Measuring urine osmolality and electrolyte-free water clearance hour by hour separates the two contributions. |
| P4 | A proactive DDAVP clamp does not so much lower the mean correction rate as **collapse the variance between phenotypes.** | Compare the *variance* of the 24-hour rise between a clamped group and a control group. The variance, not the mean, should be the primary endpoint. |
| P5 | Urea does not reduce the osmotic injury itself (Ω); it acts on the BBB and microglial arms only. | If, in an animal experiment with the correction rate strictly matched, urea reduces an *early* marker of astrocyte death, this model is wrong. |
| P6 | Relowering has a deadline, and it is set by the time constant of astrocyte death. | A cohort stratified by the time of relowering after overcorrection. |
| P7 | The risk contributed by hypokalaemia decomposes into the Edelman term (raising the sodium) and the Na⁺-coupled transport term (slowing osmolyte reaccumulation), and the latter is the larger. | Compare with and without potassium repletion while holding the rise in sodium matched. |
| P8 | In acute hyponatraemia (<48 h), Ω never becomes positive at any correction rate. | A case of ODS arising in the acute phenotype is a counterexample to this prediction (see entry 108 — the mechanism there may in fact be a different one). |

---

## Appendix B — How the citations were verified

Every PMID was confirmed by the following procedure.

```
esearch.fcgi?db=pubmed&term=<author[au] AND title-word[ti] AND year[dp]>
  -> esummary.fcgi?db=pubmed&id=<pmid>
  -> compare the first author / year / journal / title returned against the citation in the text
```

Where the search failed to return the intended paper (entries 35 and 44, for example) no PMID was attached and
the entry was left pointing at the related series instead. The principle followed was **never to guess at a PMID
that has not been confirmed**.
