# Normal pressure hydrocephalus (iNPH) QSP model — references

> **Verification method.** Every PMID in sections 1–18 below was looked up through the NCBI PubMed
> E-utilities API (`esearch` + `esummary`) and **its authors, year, journal and title confirmed to
> match the actual record**. No citation here was written from memory.
> Section 19 does the opposite: it states **which parts of the model are not supported by the literature** —
> and the parameters that the sensitivity analysis showed actually own the answer
> (`kW_out`, `Hcol_cm`, `kAQ_rec`) are precisely the ones on that list.
>
> Format: `author year journal — title [PMID]`

---

## 1. Disease concept · diagnostic criteria · guidelines

1. Adams RD 1965 *N Engl J Med* — Symptomatic occult hydrocephalus with "normal" cerebrospinal-fluid pressure: a treatable syndrome. [PMID 14303656](https://pubmed.ncbi.nlm.nih.gov/14303656/) — The original description. The very paradox that is this model's starting point: "the pressure is normal yet the disease is there".
2. Relkin N 2005 *Neurosurgery* — Diagnosing idiopathic normal-pressure hydrocephalus. [PMID 16160425](https://pubmed.ncbi.nlm.nih.gov/16160425/)
3. Marmarou A 2005 *Neurosurgery* — The value of supplemental prognostic tests for the preoperative assessment of idiopathic normal-pressure hydrocephalus. [PMID 16160426](https://pubmed.ncbi.nlm.nih.gov/16160426/)
4. Marmarou A 2005 *Acta Neurochir Suppl* — Guidelines for management of idiopathic normal pressure hydrocephalus: progress to date. [PMID 16463856](https://pubmed.ncbi.nlm.nih.gov/16463856/)
5. Klinge P 2005 *Neurosurgery* — Outcome of shunting in idiopathic normal-pressure hydrocephalus and the value of outcome assessment in shunted patients. [PMID 16160428](https://pubmed.ncbi.nlm.nih.gov/16160428/)
6. Ishikawa M 2008 *Neurol Med Chir (Tokyo)* — Guidelines for management of idiopathic normal pressure hydrocephalus. [PMID 18408356](https://pubmed.ncbi.nlm.nih.gov/18408356/)
7. Mori E 2012 *Neurol Med Chir (Tokyo)* — Guidelines for management of idiopathic normal pressure hydrocephalus: second edition. [PMID 23183074](https://pubmed.ncbi.nlm.nih.gov/23183074/)
8. Nakajima M 2021 *Neurol Med Chir (Tokyo)* — Guidelines for management of idiopathic normal pressure hydrocephalus (third edition). [PMID 33455998](https://pubmed.ncbi.nlm.nih.gov/33455998/) — The primary source for the diagnostic criteria and the investigation algorithm.
9. Tsakanikas D 2007 *Semin Neurol* — Normal pressure hydrocephalus. [PMID 17226742](https://pubmed.ncbi.nlm.nih.gov/17226742/)
10. Rovira À 2026 *Can Assoc Radiol J* — Idiopathic normal pressure hydrocephalus: a comprehensive review. [PMID 42239999](https://pubmed.ncbi.nlm.nih.gov/42239999/)
11. Schmidt E 2025 *J Neurosurg Sci* — Treatment of iNPH: novel insights. [PMID 40045806](https://pubmed.ncbi.nlm.nih.gov/40045806/)
12. Taddei G 2026 *Eur J Health Econ* — Economic analysis and healthcare implications of underdiagnosed idiopathic normal pressure hydrocephalus. [PMID 42201617](https://pubmed.ncbi.nlm.nih.gov/42201617/)

## 2. Quantitative theory of CSF dynamics — the Davson equation, the Marmarou pressure-volume curve, outflow resistance

13. Marmarou A 1975 *J Neurosurg* — Compartmental analysis of compliance and outflow resistance of the cerebrospinal fluid system. [PMID 1181384](https://pubmed.ncbi.nlm.nih.gov/1181384/) — **The model's entire hydraulic block follows the form of this paper**: the exponential pressure-volume curve `C = 1/(E1·(P−P0))` and the outflow resistance `Rout`.
14. Yamada S 2023 *Neurol Med Chir (Tokyo)* — Cerebrospinal fluid production and absorption and ventricular enlargement mechanisms in hydrocephalus. [PMID 36858632](https://pubmed.ncbi.nlm.nih.gov/36858632/) — `If0 = 0.35 mL/min`, and the distribution of absorption routes.
15. Hansson W 2026 *Neurosurgery* — Intracranial pressure dynamics and cerebrospinal fluid outflow resistance in the elderly. [PMID 41854330](https://pubmed.ncbi.nlm.nih.gov/41854330/) — The `Rout` distribution in a normal elderly group. The basis for `Rout_norm = 9`.
16. Qvarlander S 2014 *Med Biol Eng Comput* — CSF dynamic analysis of a predictive pulsatility-based infusion test for normal pressure hydrocephalus. [PMID 24151060](https://pubmed.ncbi.nlm.nih.gov/24151060/) — How to obtain `Rout` and pulsatility together from an infusion test.
17. Lalou AD 2020 *Fluids Barriers CNS* — Cerebrospinal fluid dynamics in non-acute post-traumatic ventriculomegaly. [PMID 32228689](https://pubmed.ncbi.nlm.nih.gov/32228689/) — Clinical interpretation of `Rout`, RAP and compensatory reserve.
18. Zeiler FA 2019 *Acta Neurochir (Wien)* — Compensatory-reserve-weighted intracranial pressure versus intracranial pressure for outcome association. [PMID 31053909](https://pubmed.ncbi.nlm.nih.gov/31053909/) — The argument that compensatory reserve carries more information than the mean pressure.
19. Lindstrøm EK 2018 *Neuroimage Clin* — Cerebrospinal fluid volumetric net flow rate and direction in idiopathic normal pressure hydrocephalus. [PMID 30238917](https://pubmed.ncbi.nlm.nih.gov/30238917/)
20. Rosenberg GA 1982 *Brain Res* — The effect of increased CSF pressure on interstitial fluid flow during ventriculocisternal perfusion. [PMID 7055690](https://pubmed.ncbi.nlm.nih.gov/7055690/) — Rising pressure → interstitial fluid flow; the experimental basis for the `Wpv` compartment.
21. James AE Jr 1975 *J Neurol Sci* — Pathophysiology of chronic communicating hydrocephalus in dogs. [PMID 1172782](https://pubmed.ncbi.nlm.nih.gov/1172782/) — The classic experiment on transependymal outflow and ventricular enlargement.

## 3. The pulsatility hypothesis — the claim that the waveform, not the mean pressure, is the cause

22. Eide PK 2006 *Acta Neurochir (Wien)* — Intracranial pulse pressure amplitude levels determined during preoperative assessment of subjects with possible idiopathic normal pressure hydrocephalus. [PMID 17039303](https://pubmed.ncbi.nlm.nih.gov/17039303/) — **`E1_init` and `Vp_init` were calibrated to the AMP distribution of this paper (normal ~2, iNPH >4 mmHg).**
23. Eide PK 2010 *Cerebrospinal Fluid Res* — Cerebrospinal fluid pulse pressure amplitude during lumbar infusion in idiopathic normal pressure hydrocephalus. [PMID 20205911](https://pubmed.ncbi.nlm.nih.gov/20205911/)
24. Eide PK 2008 *Med Eng Phys* — Comparison of simultaneous continuous intracranial pressure signals from sensors placed within brain parenchyma and epidural space. [PMID 17336574](https://pubmed.ncbi.nlm.nih.gov/17336574/) — The site-dependence of the AMP measurement itself. Why this model reports AMP in the "supine" posture only.
25. Eide PK 2018 *J Neurosurg* — The pathophysiology of chronic noncommunicating hydrocephalus: lessons from continuous intracranial pressure monitoring. [PMID 28799879](https://pubmed.ncbi.nlm.nih.gov/28799879/)
26. Park EH 2012 *J Neurosurg* — Impaired pulsation absorber mechanism in idiopathic normal pressure hydrocephalus. [PMID 23061391](https://pubmed.ncbi.nlm.nih.gov/23061391/) — Failure of intracranial pulsation absorption. The mechanism `kappa_tm` sets out to represent.
27. Bateman GA 2002 *Neuroradiology* — Pulse-wave encephalopathy: a comparative study of the hydrodynamics of leukoaraiosis and normal-pressure hydrocephalus. [PMID 12221445](https://pubmed.ncbi.nlm.nih.gov/12221445/)
28. Bateman GA 2003 *Neuroradiology* — The reversibility of reduced cortical vein compliance in normal-pressure hydrocephalus following shunt insertion. [PMID 12592485](https://pubmed.ncbi.nlm.nih.gov/12592485/) — **The observation that compliance recovers over several months after shunting. The basis for `kE1_rec`.**
29. Bateman GA 2004 *Med Hypotheses* — Pulse wave encephalopathy: a spectrum hypothesis incorporating Alzheimer's disease, vascular dementia and normal pressure hydrocephalus. [PMID 14962623](https://pubmed.ncbi.nlm.nih.gov/14962623/)
30. Bateman GA 2016 *Fluids Barriers CNS* — A comparison between the pathophysiology of multiple sclerosis and normal pressure hydrocephalus. [PMID 27658732](https://pubmed.ncbi.nlm.nih.gov/27658732/)
31. Murambi RT 2025 *Fluids Barriers CNS* — Deformation of brain in normal pressure hydrocephalus is more readily associated with slow vasomotion than with the cardiac pulse. [PMID 40533775](https://pubmed.ncbi.nlm.nih.gov/40533775/) — **A reference that partly conflicts with the model**: this model hangs `kappa_tm` on the cardiac pulsation component. See section 19.

## 4. Glymphatic clearance · AQP4

32. Hasan-Olive MM 2019 *Glia* — Loss of perivascular aquaporin-4 in idiopathic normal pressure hydrocephalus. [PMID 30306658](https://pubmed.ncbi.nlm.nih.gov/30306658/) — **The direct basis for `AQ_init = 0.45` (normal 0.85).**
33. Eide PK 2020 *Brain Commun* — MRI biomarkers of cerebrospinal fluid tracer dynamics in idiopathic normal pressure hydrocephalus. [PMID 33381757](https://pubmed.ncbi.nlm.nih.gov/33381757/)
34. Jacobsen HH 2020 *Invest Ophthalmol Vis Sci* — In vivo evidence for impaired glymphatic function in the visual pathway of patients with normal pressure hydrocephalus. [PMID 33201186](https://pubmed.ncbi.nlm.nih.gov/33201186/)
35. Eide PK 2019 *Gerontol Geriatr Med* — In vivo imaging of molecular clearance from human entorhinal cortex. [PMID 31819895](https://pubmed.ncbi.nlm.nih.gov/31819895/)
36. Eide PK 2022 *J Cereb Blood Flow Metab* — Altered glymphatic enhancement of cerebrospinal fluid tracer in individuals with chronic poor sleep quality. [PMID 35350917](https://pubmed.ncbi.nlm.nih.gov/35350917/) — **The human evidence for the sleep → glymphatic coupling (`sleep_base`, the melatonin scenario S19).**
37. Yang Y 2025 *Eur J Neurol* — Alterations of glymphatic system before and after shunt surgery in patients with idiopathic normal pressure hydrocephalus. [PMID 40365713](https://pubmed.ncbi.nlm.nih.gov/40365713/) — Glymphatic recovery after shunting. The test target for the prediction of section 9.
38. Mossige I 2026 *Radiology* — Comparing glymphatic function measures: DTI-ALPS and CSF tracer dynamics. [PMID 41631990](https://pubmed.ncbi.nlm.nih.gov/41631990/)
39. Broggi M 2025 *Neurol Sci* — Implications of the glymphatic system in the diagnostic and surgical workup of normal pressure hydrocephalus. [PMID 40524080](https://pubmed.ncbi.nlm.nih.gov/40524080/)
40. He W 2025 *Neurosurg Rev* — Association between perivascular spaces, DTI-derived indices, and choroid plexus with ventriculomegaly. [PMID 41081971](https://pubmed.ncbi.nlm.nih.gov/41081971/)
41. Zahran A 2026 *CNS Neurosci Ther* — Glymphatic system dysfunction in central nervous system diseases. [PMID 41792880](https://pubmed.ncbi.nlm.nih.gov/41792880/)
42. Herz J 2018 *Methods Mol Biol* — Morphological and functional analysis of CNS-associated lymphatics. [PMID 30242757](https://pubmed.ncbi.nlm.nih.gov/30242757/) — The meningeal lymphatic outflow route.

## 5. Imaging morphology — DESH, Evans index, callosal angle

43. Yamada S 2026 *Fluids Barriers CNS* — CT-based automatic segmentation of key CSF regions for detecting DESH. [PMID 42337617](https://pubmed.ncbi.nlm.nih.gov/42337617/)
44. Skalický P 2021 *J Clin Neurosci* — Role of DESH, callosal angle and cingulate sulcus sign in prediction of gait responsiveness after shunting. [PMID 33334664](https://pubmed.ncbi.nlm.nih.gov/33334664/) — The clinical meaning of the `DESH` multiplier.
45. Rohatgi S 2024 *J Comput Assist Tomogr* — Correlating Evans index, callosal angle, and lateral ventricle volume with gait response outcomes. [PMID 38657140](https://pubmed.ncbi.nlm.nih.gov/38657140/)
46. Rohatgi S 2024 *Cureus* — Predicting gait speed improvement in idiopathic normal pressure hydrocephalus patients: the role of Evans index. [PMID 38975434](https://pubmed.ncbi.nlm.nih.gov/38975434/)
47. Selcuk M 2026 *Surg Radiol Anat* — Normative Evans index values in the adult population: an age- and sex-stratified MRI study. [PMID 42262516](https://pubmed.ncbi.nlm.nih.gov/42262516/) — The normal baseline for `evans_index()`.
48. Neikter J 2020 *AJNR Am J Neuroradiol* — Ventricular volume is more strongly associated with clinical improvement than the Evans index after shunting. [PMID 32527841](https://pubmed.ncbi.nlm.nih.gov/32527841/) — **Why the model does not use the Evans index as a prognostic marker.**
49. Virhammar J 2019 *J Neurosurg* — Increase in callosal angle and decrease in ventricular volume after shunt surgery. [PMID 29393749](https://pubmed.ncbi.nlm.nih.gov/29393749/) — The magnitude of the morphological change after shunting. The calibration target for `f_plastic` (the plastic expansion fraction).
50. Holmgren RT 2026 *Fluids Barriers CNS* — Ventricular volumetry in relation to clinical response and overdrainage after shunt surgery. [PMID 42381066](https://pubmed.ncbi.nlm.nih.gov/42381066/)
51. Yeh PY 2026 *J Neuroradiol* — Predicting ventricular volume reduction and overdrainage in idiopathic normal pressure hydrocephalus. [PMID 41213359](https://pubmed.ncbi.nlm.nih.gov/41213359/)
52. Sahuquillo J 2026 *Biomedicines* — Redefining idiopathic normal pressure hydrocephalus using AI-driven brain volumetry. [PMID 41898322](https://pubmed.ncbi.nlm.nih.gov/41898322/)
53. Barough SS 2025 *medRxiv* — Disproportionately elevated sulcal index (DESI). [PMID 41404281](https://pubmed.ncbi.nlm.nih.gov/41404281/) — *Preprint (not yet peer-reviewed)*.

## 6. Cerebral blood flow · periventricular white matter

54. Virhammar J 2017 *AJNR Am J Neuroradiol* — Arterial spin-labeling perfusion MR imaging demonstrates regional CBF decrease in idiopathic normal pressure hydrocephalus. [PMID 28860216](https://pubmed.ncbi.nlm.nih.gov/28860216/) — **The calibration basis for `CBF0` and for the magnitude of the periventricular perfusion deficit in iNPH.**
55. Ziegelitz D 2014 *J Magn Reson Imaging* — Cerebral perfusion measured by DSC MRI is reduced in patients with idiopathic normal pressure hydrocephalus. [PMID 24006249](https://pubmed.ncbi.nlm.nih.gov/24006249/)
56. Ziegelitz D 2016 *J Cereb Blood Flow Metab* — Pre- and postoperative cerebral blood flow changes in patients with idiopathic normal pressure hydrocephalus. [PMID 26661191](https://pubmed.ncbi.nlm.nih.gov/26661191/) — Recovery of perfusion after shunting.
57. Tanaka A 1997 *Neurosurgery* — Cerebral blood flow and autoregulation in normal pressure hydrocephalus. [PMID 9179888](https://pubmed.ncbi.nlm.nih.gov/9179888/) — The basis for the `Autoreg` state variable.
58. Tuniz F 2017 *Fluids Barriers CNS* — The role of perfusion and diffusion MRI in the assessment of patients affected by probable idiopathic normal pressure hydrocephalus. [PMID 28899431](https://pubmed.ncbi.nlm.nih.gov/28899431/)
59. Ades-Aron B 2018 *AJNR Am J Neuroradiol* — Diffusional kurtosis along the corticospinal tract in adult normal pressure hydrocephalus. [PMID 30385473](https://pubmed.ncbi.nlm.nih.gov/30385473/) — **The anatomical basis for the leg corticospinal tract being the pathway closest to the ventricle (`CST_leg`).**
60. Kanno S 2017 *Fluids Barriers CNS* — A change in brain white matter after shunt surgery in idiopathic normal pressure hydrocephalus: a tract-based spatial statistics study. [PMID 28132644](https://pubmed.ncbi.nlm.nih.gov/28132644/) — The existence of a recoverable component.
61. Di Curzio DL 2013 *Exp Neurol* — Reduced subventricular zone proliferation and white matter damage in juvenile ferrets with kaolin-induced hydrocephalus. [PMID 23769908](https://pubmed.ncbi.nlm.nih.gov/23769908/)
62. Zadka Y 2023 *J Appl Physiol* — Mechanisms of reduced cerebral blood flow in cerebral edema and elevated intracranial pressure. [PMID 36603049](https://pubmed.ncbi.nlm.nih.gov/36603049/) — The mechanistic basis for `kCBF_W` (interstitial oedema → hypoperfusion).

## 7. CSF biomarkers — and the dilution artefact

63. Jingami N 2015 *J Alzheimers Dis* — Idiopathic normal pressure hydrocephalus has a different cerebrospinal fluid biomarker profile from Alzheimer's disease. [PMID 25428256](https://pubmed.ncbi.nlm.nih.gov/25428256/) — **The key observation of section 9: in iNPH, Aβ42 and p-tau are *both* low.**
64. Jingami N 2019 *J Alzheimers Dis* — Two-point dynamic observation of Alzheimer's disease cerebrospinal fluid biomarkers in idiopathic normal pressure hydrocephalus. [PMID 31561378](https://pubmed.ncbi.nlm.nih.gov/31561378/)
65. Miyajima M 2013 *PLoS One* — Leucine-rich α2-glycoprotein is a novel biomarker of neurodegenerative disease in human cerebrospinal fluid. [PMID 24058569](https://pubmed.ncbi.nlm.nih.gov/24058569/) — The `A_lrg` compartment.
66. Grønning R 2023 *Fluids Barriers CNS* — Association between ventricular CSF biomarkers and outcome after shunt surgery in idiopathic normal pressure hydrocephalus. [PMID 37880775](https://pubmed.ncbi.nlm.nih.gov/37880775/) — **The value differs with the sampling site (ventricular versus lumbar) — the limitation of the model using a single CSF compartment.**
67. Migliorati K 2021 *Neurol Res* — P-tau as prognostic marker in long-term follow-up for patients with shunted iNPH. [PMID 33059546](https://pubmed.ncbi.nlm.nih.gov/33059546/)
68. Kanemoto H 2023 *Int Psychogeriatr* — Cerebrospinal fluid amyloid beta and response of cognition to a tap test in idiopathic normal pressure hydrocephalus. [PMID 34399871](https://pubmed.ncbi.nlm.nih.gov/34399871/)
69. Cihlo M 2025 *Neurosurg Rev* — Value of biomarkers in the prediction of shunt responsivity in patients with normal pressure hydrocephalus. [PMID 40457021](https://pubmed.ncbi.nlm.nih.gov/40457021/)
70. Kwiecień A 2026 *Molecules* — Decoding the CSF proteomic signature of idiopathic normal pressure hydrocephalus: a systematic review. [PMID 42451686](https://pubmed.ncbi.nlm.nih.gov/42451686/)
71. Lolansen SD 2021 *Dis Markers* — Inflammatory markers in cerebrospinal fluid from patients with hydrocephalus: a systematic literature review. [PMID 33613789](https://pubmed.ncbi.nlm.nih.gov/33613789/)

## 8. Concomitant Alzheimer pathology — what a responder and a non-responder really are

72. Ye BS 2025 *Alzheimers Dement* — Degenerative pathologies on cortical biopsy, dopaminergic depletion, and shunt efficacy in iNPH. [PMID 41366839](https://pubmed.ncbi.nlm.nih.gov/41366839/) — **The basis for the `APOE`/`Ab_plq` terms being the model's only non-responder mechanism.**
73. Luikku AJ 2024 *J Neuropathol Exp Neurol* — Deep learning assisted quantitative analysis of Aβ and microglia in patients with idiopathic normal pressure hydrocephalus. [PMID 39101555](https://pubmed.ncbi.nlm.nih.gov/39101555/)
74. Greenberg ABW 2024 *Cereb Cortex* — Utility of cortical tissue analysis in normal pressure hydrocephalus. [PMID 38275188](https://pubmed.ncbi.nlm.nih.gov/38275188/)
75. Mattoli MV 2020 *Int J Mol Sci* — Usefulness of brain PET with different tracers in the evaluation of patients with idiopathic normal pressure hydrocephalus. [PMID 32906629](https://pubmed.ncbi.nlm.nih.gov/32906629/)

## 9. Diagnostic perturbations — tap test, continuous lumbar drainage, infusion test

76. Rydja J 2021 *Fluids Barriers CNS* — Evaluating the cerebrospinal fluid tap test with the Hellström iNPH scale. [PMID 33827613](https://pubmed.ncbi.nlm.nih.gov/33827613/)
77. Liu C 2021 *Clin Neurol Neurosurg* — A pilot study of multiple time points and multidomain assessment in cerebrospinal fluid tap test. [PMID 34749022](https://pubmed.ncbi.nlm.nih.gov/34749022/) — **The timing of assessment changes the result — directly matching the time-constant mismatch argument of section 6.**
78. Mládek A 2022 *Neurosurgery* — Prediction of shunt responsiveness in suspected patients with normal pressure hydrocephalus using the lumbar infusion test. [PMID 35080523](https://pubmed.ncbi.nlm.nih.gov/35080523/)
79. Akar K 2026 *Fluids Barriers CNS* — Sensitivity of physiotherapy-based clinical tests in detecting change in gait and balance performance. [PMID 41703573](https://pubmed.ncbi.nlm.nih.gov/41703573/) — **The evidence that test sensitivity depends on the observer's decision threshold (section 6).**
80. Wolfsegger T 2017 *J Neurol Sci* — Cognitive impairment predicts worse short-term response to spinal tap test in normal pressure hydrocephalus. [PMID 28716246](https://pubmed.ncbi.nlm.nih.gov/28716246/)
81. Kudelić N 2023 *Front Neurol* — Predictive value of spinal CSF volume in the preoperative assessment of patients with idiopathic normal-pressure hydrocephalus. [PMID 37869132](https://pubmed.ncbi.nlm.nih.gov/37869132/)
82. Cai H 2025 *J Med Internet Res* — Predictive value of digital neuropsychological and gait assessments on shunt outcome. [PMID 41289581](https://pubmed.ncbi.nlm.nih.gov/41289581/)
83. Guarracino I 2025 *Brain Sci* — Real-time neuropsychological testing during infusion in hydrocephalus. [PMID 39851404](https://pubmed.ncbi.nlm.nih.gov/39851404/)

## 10. Randomised shunt trials — valve pressure and the operation itself

84. Boon AJ 1998 *J Neurosurg* — Dutch Normal-Pressure Hydrocephalus Study: randomized comparison of low- and medium-pressure shunts. [PMID 9488303](https://pubmed.ncbi.nlm.nih.gov/9488303/) — **The trial the titration map of section 5 sets out to reproduce. The double edge that a lower opening pressure is more effective but overdrainage complications increase.**
85. Lemcke J 2013 *J Neurol Neurosurg Psychiatry* — Safety and efficacy of gravitational shunt valves in patients with idiopathic normal pressure hydrocephalus (SVASONA). [PMID 23457222](https://pubmed.ncbi.nlm.nih.gov/23457222/) — **The randomised trial of the gravitational assist device. The key evidence for sections 3 and 5.**
86. Lemcke J 2012 *Acta Neurochir Suppl* — On the method of a randomised comparison of programmable valves with and without gravitational units. [PMID 22327702](https://pubmed.ncbi.nlm.nih.gov/22327702/)
87. Lemcke J 2010 *Acta Neurochir Suppl* — Is it possible to minimize overdrainage complications with gravitational units in patients with idiopathic normal pressure hydrocephalus? [PMID 19812931](https://pubmed.ncbi.nlm.nih.gov/19812931/)
88. Luciano M 2023 *Neurosurgery* — Placebo-controlled effectiveness of idiopathic normal pressure hydrocephalus shunting: a randomized pilot. [PMID 36700738](https://pubmed.ncbi.nlm.nih.gov/36700738/)
89. Luciano MG 2025 *N Engl J Med* — A randomized trial of shunting for idiopathic normal-pressure hydrocephalus. [PMID 40960253](https://pubmed.ncbi.nlm.nih.gov/40960253/) — **The placebo (closed-valve) controlled trial. The reference point for the section 12 self-criticism that the model overestimates the shunt effect.**
90. Saper CB 2026 *Ann Neurol* — Time to reconsider the value of shunting procedures for "idiopathic normal pressure hydrocephalus". [PMID 41741941](https://pubmed.ncbi.nlm.nih.gov/41741941/) — **A commentary in direct tension with the model's conclusion. It must be read alongside.**
91. Klinge P 2012 *Acta Neurol Scand* — One-year outcome in the European multicentre study on iNPH. [PMID 22571428](https://pubmed.ncbi.nlm.nih.gov/22571428/) — The actual response rate, 60–80%.
92. Salih A 2024 *EClinicalMedicine* — The effectiveness of various CSF diversion surgeries in idiopathic normal pressure hydrocephalus: a systematic review. [PMID 39539993](https://pubmed.ncbi.nlm.nih.gov/39539993/)
93. Brenner LBO 2026 *Neurosurgery* — The Lumboperitoneal Shunt Study: a systematic review and single-arm meta-analysis of 2696 patients. [PMID 42294933](https://pubmed.ncbi.nlm.nih.gov/42294933/)
94. Burrows EJ 2026 *Neurosurg Pract* — Clinical subgroups and treatment outcomes in idiopathic normal pressure hydrocephalus. [PMID 42282951](https://pubmed.ncbi.nlm.nih.gov/42282951/)

## 11. Shunt hydraulics · hydrostatics · overdrainage — section 3 of this model

95. de Jong DA 2000 *Acta Neurochir (Wien)* — Hydrostatic and hydrodynamic considerations in shunted normal pressure hydrocephalus. [PMID 10819253](https://pubmed.ncbi.nlm.nih.gov/10819253/) — **The reference the hydrostatic column calculation of section 3 follows. The basis for `Hcol_cm`.**
96. Kajimoto Y 2000 *J Neurosurg* — Posture-related changes in the pressure environment of the ventriculoperitoneal shunt system. [PMID 11014539](https://pubmed.ncbi.nlm.nih.gov/11014539/) — **Direct measurement showing that the driving pressure rises by tens of cmH₂O on standing. The justification for the `f_up` weighting scheme.**
97. Medow JE 2012 *J Neurosurg Pediatr* — Posture-independent piston valve: a novel valve mechanism that actuates based on intracranial pressure alone. [PMID 22208323](https://pubmed.ncbi.nlm.nih.gov/22208323/)
98. Corin AS 2026 *Neurochirurgie* — Effectiveness of anti-siphon devices in CSF shunts for preventing overdrainage in normal pressure hydrocephalus. [PMID 42486335](https://pubmed.ncbi.nlm.nih.gov/42486335/) — `asd_eff`.
99. Hafizka Y 2026 *World Neurosurg* — Effect of anti-siphon devices on postoperative outcomes in idiopathic normal pressure hydrocephalus. [PMID 41941959](https://pubmed.ncbi.nlm.nih.gov/41941959/)
100. Meier U 2013 *Neurosurgery* — Predictors of subsequent overdrainage and clinical outcomes after ventriculoperitoneal shunting for idiopathic normal pressure hydrocephalus. [PMID 24257332](https://pubmed.ncbi.nlm.nih.gov/24257332/) — Atrophy (`atrophy`) as a risk factor for overdrainage.
101. Chung K 2026 *JAAPA* — Ventriculoperitoneal shunt management in patients with normal pressure hydrocephalus and cerebral atrophy. [PMID 41662144](https://pubmed.ncbi.nlm.nih.gov/41662144/)
102. Kawahara T 2022 *Surg Neurol Int* — Dural sac shrinkage signs on spinal MRI indicate overdrainage after lumboperitoneal shunt. [PMID 35855156](https://pubmed.ncbi.nlm.nih.gov/35855156/)
103. Younes B 2026 *Childs Nerv Syst* — Incidence and management of overdrainage in hydrocephalus patients treated with adjustable valves. [PMID 41758240](https://pubmed.ncbi.nlm.nih.gov/41758240/)
104. Chelmis F 2026 *Cureus* — Diagnostic challenges and surgical outcomes of Miyazaki syndrome. [PMID 42460193](https://pubmed.ncbi.nlm.nih.gov/42460193/) — The extreme phenotype of overdrainage.
105. Iimori T 2026 *Clin Neurol Neurosurg* — Safety of continuing antithrombotic therapy during lumboperitoneal shunting for idiopathic normal pressure hydrocephalus. [PMID 42068899](https://pubmed.ncbi.nlm.nih.gov/42068899/) — The `antithrombotic` risk coefficient.

## 12. Valve choice · titration · revision

106. Reis RC 2026 *J Neurol Surg A Cent Eur Neurosurg* — Optimizing initial shunt pressure in idiopathic normal pressure hydrocephalus. [PMID 41571242](https://pubmed.ncbi.nlm.nih.gov/41571242/) — **The clinical counterpart of the titration map of section 5.**
107. Colonna S 2026 *World Neurosurg* — Valve selection and long-term outcomes in idiopathic normal pressure hydrocephalus. [PMID 41720264](https://pubmed.ncbi.nlm.nih.gov/41720264/)
108. Ahmed M 2023 *Clin Neurol Neurosurg* — Fixed versus adjustable differential pressure valves in idiopathic normal pressure hydrocephalus. [PMID 37209623](https://pubmed.ncbi.nlm.nih.gov/37209623/)
109. Katiyar V 2021 *Neurol India* — Comparison of programmable and non-programmable shunts for normal pressure hydrocephalus: a meta-analysis. [PMID 35102997](https://pubmed.ncbi.nlm.nih.gov/35102997/)
110. Oksa S 2026 *Acta Neurochir (Wien)* — Reduced risk of shunt revision with adjustable valves: a population-based cohort study over three decades. [PMID 41634436](https://pubmed.ncbi.nlm.nih.gov/41634436/) — The revision rate (`k_occl`, scenario S18).

## 13. Endoscopic third ventriculostomy (ETV)

111. Scalia G 2025 *Brain Sci* — The effects of endoscopic third ventriculostomy versus ventriculoperitoneal shunt on neuropsychological outcomes. [PMID 40426679](https://pubmed.ncbi.nlm.nih.gov/40426679/)
112. Harbaugh TD 2025 *Cureus* — Ventriculoperitoneal shunting versus endoscopic third ventriculostomy for the surgical management of idiopathic normal pressure hydrocephalus. [PMID 40034632](https://pubmed.ncbi.nlm.nih.gov/40034632/) — **The clinical comparison for S20 (the prediction that ETV is ineffective).**
113. Mohamed B 2024 *Cureus* — A single-centre experience of the management and surgical outcomes of late-onset idiopathic aqueductal stenosis. [PMID 38868257](https://pubmed.ncbi.nlm.nih.gov/38868257/) — The contrast that ETV does work in obstructive hydrocephalus.

## 14. Drug treatment — the PK/PD of acetazolamide

114. Kunka RL 1979 *J Pharm Sci* — Nonlinear model for acetazolamide. [PMID 423125](https://pubmed.ncbi.nlm.nih.gov/423125/) — **The direct basis for the saturable erythrocyte binding structure (`az_Bmax`, `az_kon`, `az_koff`).**
115. Maren TH 1979 *J Appl Physiol* — Effect of varying CO₂ equilibria on rates of HCO₃⁻ formation in cerebrospinal fluid. [PMID 118142](https://pubmed.ncbi.nlm.nih.gov/118142/) — The coupling of carbonic anhydrase to CSF formation.
116. Vogh BP 1975 *Am J Physiol* — Sodium, chloride, and bicarbonate movement from plasma to cerebrospinal fluid. [PMID 803792](https://pubmed.ncbi.nlm.nih.gov/803792/) — The animal evidence for `az_Emax = 0.50`.
117. Alperin N 2014 *Neurology* — Low-dose acetazolamide reverses periventricular white matter hyperintensities in iNPH. [PMID 24634454](https://pubmed.ncbi.nlm.nih.gov/24634454/)
118. Gilbert GJ 2014 *Neurology* — Low-dose acetazolamide reverses periventricular white matter hyperintensities in iNPH (correspondence). [PMID 25367060](https://pubmed.ncbi.nlm.nih.gov/25367060/)
119. Virhammar J 2026 *Lancet Neurol* — Safety, tolerability, and efficacy of acetazolamide in idiopathic normal pressure hydrocephalus (DRAIN). [PMID 42127932](https://pubmed.ncbi.nlm.nih.gov/42127932/) — **The clinical test target for the conclusion of section 8 (the chronic effect is small and the adverse effects are dose-dependent). The result of this randomised trial decides the interpretation of the model's S12/S13.**

## 15. The physiology of CSF secretion — at the molecular level

120. Jensen DB 2025 *Adv Sci* — The Na⁺,K⁺,2Cl⁻ cotransporter, not aquaporin 1, sustains cerebrospinal fluid secretion. [PMID 39692709](https://pubmed.ncbi.nlm.nih.gov/39692709/) — **The basis for the `NKCC1_cp` node and for the loop diuretic (`bu_Emax`) target. Why the map places NKCC1 above AQP1.**
121. Deffner F 2022 *Cell Mol Life Sci* — Aquaporin-4 expression in the human choroid plexus. [PMID 35072772](https://pubmed.ncbi.nlm.nih.gov/35072772/)

## 16. Clinical symptoms · assessment scales

122. Hellström P 2012 *Acta Neurol Scand* — A new scale for assessment of severity and outcome in iNPH. [PMID 22587624](https://pubmed.ncbi.nlm.nih.gov/22587624/)
123. Passaretti M 2023 *Mov Disord Clin Pract* — Gait analysis in idiopathic normal pressure hydrocephalus: a meta-analysis. [PMID 38026510](https://pubmed.ncbi.nlm.nih.gov/38026510/) — **The calibration basis for `G_max` and for the iNPH gait speed baseline (≈0.6 m/s).**
124. Mills R 2026 *Fluids Barriers CNS* — Three-dimensional kinematic gait signatures of idiopathic normal pressure hydrocephalus. [PMID 42174608](https://pubmed.ncbi.nlm.nih.gov/42174608/)
125. Ekblom M 2025 *Fluids Barriers CNS* — Blinded gait assessment in idiopathic normal pressure hydrocephalus: reliability and correlation with clinical scales. [PMID 40926257](https://pubmed.ncbi.nlm.nih.gov/40926257/)
126. Lander J 2026 *Fluids Barriers CNS* — Slow unsteady gait in a population-based cohort: links to ventriculomegaly and iNPH-related imaging markers. [PMID 42286688](https://pubmed.ncbi.nlm.nih.gov/42286688/)
127. Sakakibara R 2008 *Neurourol Urodyn* — Mechanism of bladder dysfunction in idiopathic normal pressure hydrocephalus. [PMID 18092331](https://pubmed.ncbi.nlm.nih.gov/18092331/) — **The evidence for the frontal–pontine micturition circuit of the `Urin` compartment.**
128. Krzastek SC 2017 *Neurourol Urodyn* — Characterization of lower urinary tract symptoms in patients with idiopathic normal pressure hydrocephalus. [PMID 27490149](https://pubmed.ncbi.nlm.nih.gov/27490149/)
129. Rovčanin B 2025 *Diagnostics* — Time course of symptoms in normal-pressure hydrocephalus: a systematic review. [PMID 40722526](https://pubmed.ncbi.nlm.nih.gov/40722526/) — The clinical background to section 7 (early versus delayed surgery).
130. Kimura T 2021 *World Neurosurg* — Preoperative predictive factors of short-term outcome in idiopathic normal pressure hydrocephalus. [PMID 33895373](https://pubmed.ncbi.nlm.nih.gov/33895373/)
131. Uchigami H 2022 *Clin Neurol Neurosurg* — Preoperative factors associated with shunt responsiveness in patients with idiopathic normal-pressure hydrocephalus. [PMID 36049404](https://pubmed.ncbi.nlm.nih.gov/36049404/) — **The observation that symptom duration predicts prognosis. The phenomenon `k_perm` (the irreversible conversion rate) sets out to represent.**
132. Mansour MA 2025 *J Neurosurg Case Lessons* — Beyond the triad: akinetic mutism in idiopathic normal pressure hydrocephalus. [PMID 40720908](https://pubmed.ncbi.nlm.nih.gov/40720908/)

## 17. Pathogenesis · genetics · epithelium

133. Sato H 2016 *PLoS One* — A segmental copy number loss of the SFMBT1 gene is a genetic risk for shunt-responsive idiopathic normal pressure hydrocephalus. [PMID 27861535](https://pubmed.ncbi.nlm.nih.gov/27861535/)
134. Tipton PW 2023 *Neurol Genet* — CWH43 variants are associated with disease risk and clinical phenotypic measures in patients with normal pressure hydrocephalus. [PMID 37476022](https://pubmed.ncbi.nlm.nih.gov/37476022/)
135. Yang HW 2023 *Proc Natl Acad Sci USA* — A role for mutations in AK9 and other genes affecting ependymal cells in idiopathic normal pressure hydrocephalus. [PMID 38100419](https://pubmed.ncbi.nlm.nih.gov/38100419/) — The `Epend_loss` node.
136. Piccinin CC 2025 *Mov Disord* — Genetic risk factors in normal pressure hydrocephalus: what we know and what is next. [PMID 40266017](https://pubmed.ncbi.nlm.nih.gov/40266017/)
137. Yang D 2021 *Proc Natl Acad Sci USA* — Increased plasmin-mediated proteolysis of L1CAM in a mouse model of idiopathic normal pressure hydrocephalus. [PMID 34380733](https://pubmed.ncbi.nlm.nih.gov/34380733/)
138. Botfield H 2013 *Brain* — Decorin prevents the development of juvenile communicating hydrocephalus. [PMID 23983032](https://pubmed.ncbi.nlm.nih.gov/23983032/) — **The causal experimental evidence for the TGF-β1 → arachnoid fibrosis → rising `Rout` pathway.**
139. Yan H 2016 *Brain Res* — Decorin alleviated chronic hydrocephalus via inhibiting TGF-β1/Smad/CTGF pathway after subarachnoid hemorrhage. [PMID 26556770](https://pubmed.ncbi.nlm.nih.gov/26556770/)
140. Tan Q 2017 *Brain Res* — Cannabinoid receptor 2 activation restricts fibrosis and alleviates hydrocephalus after intraventricular hemorrhage. [PMID 27769788](https://pubmed.ncbi.nlm.nih.gov/27769788/)
141. Feng Z 2020 *Neurosci Lett* — uPA alleviates kaolin-induced hydrocephalus. [PMID 32497735](https://pubmed.ncbi.nlm.nih.gov/32497735/)
142. Ishikawa M 2008 *Brain Nerve* — [Idiopathic normal pressure hydrocephalus — overviews and pathogenesis]. [PMID 18402067](https://pubmed.ncbi.nlm.nih.gov/18402067/) — *In Japanese*.

## 18. Sleep-disordered breathing · differential diagnosis · mathematical modelling

143. Riedel CS 2022 *Sleep* — Sleep-disordered breathing is frequently associated with idiopathic normal pressure hydrocephalus. [PMID 34739077](https://pubmed.ncbi.nlm.nih.gov/34739077/) — `OSA_node` → `Sleep_SWS` → glymphatic.
144. Román GC 2019 *Curr Neurol Neurosci Rep* — Sleep-disordered breathing and idiopathic normal-pressure hydrocephalus: recent pathophysiological advances. [PMID 31144048](https://pubmed.ncbi.nlm.nih.gov/31144048/)
145. Riedel CS 2024 *Fluids Barriers CNS* — Transient intracranial pressure elevations (B waves) associated with sleep. [PMID 39702226](https://pubmed.ncbi.nlm.nih.gov/39702226/)
146. Yun S 2025 *Neuroradiology* — Data-driven differentiation of idiopathic normal-pressure hydrocephalus and progressive supranuclear palsy. [PMID 41204957](https://pubmed.ncbi.nlm.nih.gov/41204957/)
147. Deng Z 2024 *Clin Neurol Neurosurg* — Evaluation of imaging indicators in differentiating idiopathic normal pressure hydrocephalus from Alzheimer's disease. [PMID 38823198](https://pubmed.ncbi.nlm.nih.gov/38823198/)
148. Netterwala A 2026 *BMJ Case Rep* — Recognising CADASIL in adults with NPH-like syndrome. [PMID 42236108](https://pubmed.ncbi.nlm.nih.gov/42236108/)
149. Suppa A 2025 *Mov Disord* — Neurophysiology of atypical parkinsonian syndromes: a study group position paper. [PMID 40356334](https://pubmed.ncbi.nlm.nih.gov/40356334/)
150. Dreyer LW 2024 *Fluids Barriers CNS* — Modeling CSF circulation and the glymphatic system during infusion using subject-specific intracranial pressures. [PMID 39407250](https://pubmed.ncbi.nlm.nih.gov/39407250/) — **The closest prior modelling study. The difference from this model is that this one carries drug PK/PD and shunt hardware as well.**
151. Chu KH 2024 *Brain Spine* — Mathematical modelling of cerebral haemodynamics and their effects on ICP. [PMID 38510619](https://pubmed.ncbi.nlm.nih.gov/38510619/)
152. Doron O 2021 *Fluids Barriers CNS* — Interactions of brain, blood, and CSF: a novel mathematical model of cerebral edema. [PMID 34530863](https://pubmed.ncbi.nlm.nih.gov/34530863/)
153. Ursino M 2010 *Ann Biomed Eng* — A model of cerebrovascular reactivity including the circle of Willis and cortical anastomoses. [PMID 20094916](https://pubmed.ncbi.nlm.nih.gov/20094916/) — The classic lineage of cerebrovascular compliance and autoregulation modelling.
154. Lee KJ 2015 *Biomed Eng Online* — Non-invasive detection of intracranial hypertension using a simplified intracranial hemo- and hydro-dynamic model. [PMID 26024843](https://pubmed.ncbi.nlm.nih.gov/26024843/)
155. Gadda G 2015 *Am J Physiol Heart Circ Physiol* — A new hemodynamic model for the study of cerebral venous outflow. [PMID 25398980](https://pubmed.ncbi.nlm.nih.gov/25398980/)
156. Lampe R 2014 *Comput Math Methods Med* — Mathematical modelling of cerebral blood circulation and cerebral autoregulation. [PMID 25126111](https://pubmed.ncbi.nlm.nih.gov/25126111/)

---

## 19. The parts that are **not** sourced from the literature

What is dangerous in a QSP model is not a parameter with no citation but **a parameter whose lack of a
citation is hidden**. Below are the items for which no directly measured value could be obtained from any
of the references above; it is the same list as section 12 of `inph_model_report.txt`.

| Parameter | Value | Status |
|---|---|---|
| `kW_out` | 1.300 /day | The clearance rate of periventricular interstitial water. **First in the sensitivity analysis (elasticity +2.27) — the parameter the answer actually hangs on, and it is an order-of-magnitude estimate.** References 20 and 21 support the existence of transependymal flow but not its rate. |
| `kAQ_rec` | 0.008 /day | The rate of AQP4 repolarisation after decompression. **Third in the sensitivity analysis.** Reference 32 shows the loss of polarity, but there is no human measurement of the recovery rate. |
| `kappa_tm` | 0.060 | **Never measured in a human ventricle.** The fraction of the ICP pulsation that crosses the mantle. **It was expected to be the most influential parameter in the whole model, yet came out tenth of fourteen in the sensitivity analysis (−0.049) — both the expectation and its refutation have been left in the report.** References 26 and 27 support the existence of the mechanism but not its magnitude. Reference 31 reports that vasomotion slower than the cardiac pulse explains the deformation better, which puts it **in tension with the very form of this term**. |
| `DESH` | 1.55 | The morphological multiplier applied to `kappa_tm`. Reference 44 shows that DESH is associated with prognosis, but the value of the multiplier is pure estimate. |
| `k_perm` | 1.0e-4 /day | The conversion rate from recoverable to permanent damage. **This one parameter determines the result of section 7 (early versus delayed surgery).** References 131 and 129 support the direction but not the rate. |
| `kflux_ab`, `kflux_tau`, `kdeg_csf` | — | The brain→CSF biomarker efflux constants. The **direction** of the prediction in section 9 hangs on this ratio, and the text says so. |
| `Hcol_cm` | 45 cmH₂O | The ventricle-to-peritoneum column, proportional to the patient's height. References 95 and 96 support the order of magnitude, but the individual value differs from patient to patient. It is the dominant term of section 3 and **second in the sensitivity analysis — in a model containing five drugs, the patient's height came out more important for gait at two years than any of them.** |
| `k_haz` | 1.5e-5 | The conversion from fluid collection volume to symptomatic-event hazard. **The subdural incidence (%) of section 5 is therefore an ordinal scale and not a calibrated incidence.** |
| `etv_dRout` | 0.06 | The claim that ETV barely lowers `Rout` in communicating hydrocephalus is a **structural argument** (the resistance lies downstream of the stoma), not a measurement. Its direction agrees with the clinical result of reference 112. |
| `az_tau_esc`, `az_f_esc` | 18 d, 0.78 | The loss of acetazolamide's chronic effect. The mechanism (renal acid-base compensation, upregulation of villous transporters) is widely accepted, but has never been measured at this rate in iNPH. Reference 119 (the DRAIN trial) will be the basis for deciding. |

### Where the model diverges from the actual literature

1. **Response rate.** This model makes almost every patient with a raised `Rout` a responder. The measured response rate of reference 91 is 60–80%, and references 89 and 90 are more sceptical still. → **The model overestimates the value of shunting in an unselected population.**
2. **Ventricular size.** The model reduces the Evans index by 0.01–0.02 after a successful shunt. References 48 and 49 report that the Evans index corresponds poorly to clinical improvement and that the change is small — the model still couples size too tightly to pressure.
3. **Acetazolamide.** Section 8 gives an acute effect of 39–43% and a chronic effect of 12%. References 117 and 118 reported improvement in white matter lesions but their samples are small, and reference 119 (the randomised trial) is the final arbiter. S12/S13 should be read as an **upper bound**.
4. **The frequency band of pulsatility.** Reference 31 reports that brain deformation is better associated with vasomotion slower than the cardiac pulse. This model hangs `kappa_tm` on `AMP` (the cardiac pulsation component). The model does not decide which of the two readings is right, and it does not hold the data with which to decide.
5. **A single CSF compartment.** Reference 66 reports that biomarkers differ between ventricular and lumbar CSF. This model holds only one CSF compartment, so it reads the tap test (lumbar) and the shunt (ventricular) at the same concentration.

---

## 20. Tools · methodology

- mrgsolve (R): <https://mrgsolve.org/>
- Graphviz: <https://graphviz.org/>
- Shiny (R): <https://shiny.posit.co/>
- NCBI E-utilities (used to verify this reference list): <https://www.ncbi.nlm.nih.gov/books/NBK25501/>

> **Disclaimer.** This model and this reference list are for educational and research purposes. They cannot be
> used for clinical decisions, prescribing, or regulatory submission. The valve titration map of section 5 in particular
> contains uncalibrated risk coefficients, and so **cannot be the basis for an actual valve
> setting.**
