# Coronary Microvascular Dysfunction (CMD / ANOCA–INOCA) — References

> Coronary Microvascular Dysfunction · Angina/Ischaemia with No Obstructive Coronary Arteries
> For the QSP model (`cmd_qsp_model_en.dot`, `cmd_mrgsolve_model.R`, `cmd_reference_model.py`):
> the basis of every parameter, structure, and calibration target.

**Citation principle.** Every PMID in this list was looked up directly with the PubMed E-utilities in this
session and its title, journal, and year confirmed. For the references used in the quantitative claims of the model (§2, §3, §11)
the full abstract was retrieved and the figures checked, and those entries also carry **the figures used
in the model**. References that could not be confirmed were not included.

---

## 1. Definitions · nomenclature · guidelines

1. Ong P, et al. **International standardization of diagnostic criteria for microvascular angina.** Int J Cardiol. 2018. [PMID 29031990](https://pubmed.ncbi.nlm.nih.gov/29031990/) — the COVADIS diagnostic criteria for microvascular angina. The endotype classification axes of the model (CFR, microvascular resistance, ACh response) follow these criteria.
2. Beltrame JF, et al. **International standardization of diagnostic criteria for vasospastic angina.** Eur Heart J. 2017. [PMID 26245334](https://pubmed.ncbi.nlm.nih.gov/26245334/) — the COVADIS criteria for vasospastic angina.
3. Vrints C, et al. **2024 ESC Guidelines for the management of chronic coronary syndromes.** Eur Heart J. 2024. [PMID 39210710](https://pubmed.ncbi.nlm.nih.gov/39210710/)
4. Vrints C, et al. **[2024 ESC guideline on chronic coronary syndromes, Italian edition].** G Ital Cardiol. 2024. [PMID 39611224](https://pubmed.ncbi.nlm.nih.gov/39611224/)
5. Gulati M, et al. **2021 AHA/ACC/ASE/CHEST/SAEM/SCCT/SCMR Guideline for the Evaluation and Diagnosis of Chest Pain.** Circulation. 2021. [PMID 34709879](https://pubmed.ncbi.nlm.nih.gov/34709879/)
6. Writing Committee Members. **2021 AHA/ACC Chest Pain Guideline (JACC).** J Am Coll Cardiol. 2021. [PMID 34756653](https://pubmed.ncbi.nlm.nih.gov/34756653/)
7. Gulati M, et al. **2021 Chest Pain Guideline: Executive Summary.** Circulation. 2021. [PMID 34709928](https://pubmed.ncbi.nlm.nih.gov/34709928/)
8. Boden WE, et al. **Myocardial Ischemic Syndromes: A New Nomenclature to Harmonize Evolving International Clinical Practice Guidelines.** Circulation. 2024. [PMID 39210827](https://pubmed.ncbi.nlm.nih.gov/39210827/)
9. Kaur G, et al. **Chest Pain in Women: Considerations From the 2021 AHA/ACC Chest Pain Guideline.** Curr Probl Cardiol. 2023. [PMID 36921653](https://pubmed.ncbi.nlm.nih.gov/36921653/)
10. Pepine CJ, et al. **ANOCA/INOCA/MINOCA: Open artery ischemia.** Am Heart J Plus. 2023. [PMID 37064505](https://pubmed.ncbi.nlm.nih.gov/37064505/)
11. Ashokprabhu ND, et al. **INOCA/ANOCA: Mechanisms and novel treatments.** Am Heart J Plus. 2023. [PMID 37377840](https://pubmed.ncbi.nlm.nih.gov/37377840/)
12. Rinaldi R, et al. **Management of angina pectoris.** Trends Cardiovasc Med. 2025. [PMID 40086653](https://pubmed.ncbi.nlm.nih.gov/40086653/)
13. Burgess S, et al. **Challenges in diagnosing coronary microvascular dysfunction and coronary vasospasm.** Cardiovasc Revasc Med. 2025. [PMID 40312200](https://pubmed.ncbi.nlm.nih.gov/40312200/)
14. Scarica V, et al. **Coronary microvascular dysfunction: pathophysiology, diagnosis, and therapeutic strategies across cardiovascular diseases.** EXCLI J. 2025. [PMID 40376434](https://pubmed.ncbi.nlm.nih.gov/40376434/)
15. Ya'Qoub L, et al. **Non-obstructive Plaque and Treatment of INOCA: More to Be Learned.** Curr Atheroscler Rep. 2022. [PMID 35781776](https://pubmed.ncbi.nlm.nih.gov/35781776/)
16. Vancheri F, et al. **Coronary Microvascular Dysfunction.** J Clin Med. 2020. [PMID 32899944](https://pubmed.ncbi.nlm.nih.gov/32899944/)
17. Smilowitz NR, et al. **Coronary Microvascular Disease in Contemporary Clinical Practice.** Circ Cardiovasc Interv. 2023. [PMID 37259860](https://pubmed.ncbi.nlm.nih.gov/37259860/)
18. Parwani P, et al. **Contemporary Diagnosis and Management of Patients with MINOCA.** Curr Cardiol Rep. 2023. [PMID 37067753](https://pubmed.ncbi.nlm.nih.gov/37067753/)
19. Boden WE, et al. **Evolving Management Paradigm for Stable Ischemic Heart Disease Patients.** J Am Coll Cardiol. 2023. [PMID 36725179](https://pubmed.ncbi.nlm.nih.gov/36725179/)

---

## 2. The backbone of this model — endotype and invasive physiology

The first reference in this section is the structural starting point of the whole model. **CFR is a ratio, and
the same low value can come from the denominator (excessive resting flow) or from the numerator (a ceiling on
hyperaemic flow).**

20. **Rahman H, et al. Coronary Microvascular Dysfunction Is Associated With Myocardial Ischemia and Abnormal Coronary Perfusion During Exercise. Circulation. 2019;140:1805–1816.** [PMID 31707835](https://pubmed.ncbi.nlm.nih.gov/31707835/)
    **Figures used in the model (calibration targets T1–T5):** 85 patients with ANOCA (78% women, 57±10 years). MVD defined as CFR<2.5 → 45 patients (53%). Of those with MVD, 62% were **functional** (hyperaemic microvascular resistance <2.5 mmHg/cm/s) and 38% **structural** (≥2.5). Resting microvascular resistance was 4.2±1.0 functional, 6.9±1.7 structural, and 7.3±2.2 mmHg/(cm/s) in controls. Inducible ischaemia (a hyperaemic subendocardial/epicardial perfusion ratio <1.0) in 82% of the MVD group against 22% of controls. Global myocardial perfusion reserve 2.01±0.41 against 2.68±0.49. Perfusion efficiency rose on exercise from 59±11 to 65±14% in controls but fell from 61±12 to 44±10% in MVD. Systolic pressure on exercise was 188±25 structural against 161±27 functional against 156±30 mmHg in controls (the basis for the model parameters `KSBP_STRUCT/FUNC/CTRL`). **Myocardial perfusion under stress and exercise perfusion efficiency were similar between the two endotypes** — a tension the model has to reproduce, and which §4 explains through the "reserve exhaustion workload".
21. Sinha A, et al. **Rethinking False Positive Exercise Electrocardiographic Stress Tests by Assessing Coronary Microvascular Function.** J Am Coll Cardiol. 2024. [PMID 38199706](https://pubmed.ncbi.nlm.nih.gov/38199706/) — 102 patients with ANOCA. Ischaemia on the exercise ECG carried **100% specificity** for CMD, and ACh flow reserve was the strongest predictor of exercise ischaemia. The premise of the model that "the exercise stress test is not a false positive; it is looking at the microvascular substrate".
22. Rahman H, et al. **Optimal Use of Vasodilators for Diagnosis of Microvascular Angina in the Cardiac Catheterization Laboratory.** Circ Cardiovasc Interv. 2020. [PMID 32519879](https://pubmed.ncbi.nlm.nih.gov/32519879/)
23. Ford TJ, et al. **Assessment of Vascular Dysfunction in Patients Without Obstructive Coronary Artery Disease: Why, How, and When.** JACC Cardiovasc Interv. 2020. [PMID 32819476](https://pubmed.ncbi.nlm.nih.gov/32819476/)
24. Suda A, et al. **Coronary Functional Abnormalities in Patients With Angina and Nonobstructive Coronary Artery Disease.** J Am Coll Cardiol. 2019. [PMID 31699275](https://pubmed.ncbi.nlm.nih.gov/31699275/)
25. Ang DTY, et al. **Phenotype-based management of coronary microvascular dysfunction.** J Nucl Cardiol. 2022. [PMID 35672569](https://pubmed.ncbi.nlm.nih.gov/35672569/)
26. Lanza GA, Crea F. **Primary coronary microvascular dysfunction: clinical presentation, pathophysiology, and management.** Circulation. 2010. [PMID 20516386](https://pubmed.ncbi.nlm.nih.gov/20516386/)
27. Crea F, et al. **Coronary microvascular dysfunction: an update.** Eur Heart J. 2014. [PMID 24366916](https://pubmed.ncbi.nlm.nih.gov/24366916/)
28. Crea F, et al. **Pathophysiology of Coronary Microvascular Dysfunction.** Circ J. 2022. [PMID 34759123](https://pubmed.ncbi.nlm.nih.gov/34759123/)
29. Del Buono MG, et al. **Coronary Microvascular Dysfunction Across the Spectrum of Cardiovascular Diseases: JACC State-of-the-Art Review.** J Am Coll Cardiol. 2021. [PMID 34556322](https://pubmed.ncbi.nlm.nih.gov/34556322/)
30. Bairey Merz CN, et al. **Treatment of coronary microvascular dysfunction.** Cardiovasc Res. 2020. [PMID 32087007](https://pubmed.ncbi.nlm.nih.gov/32087007/)
31. Nogami K, et al. **Chest pain patterns and coronary microvascular function in non-obstructive coronary artery disease.** EuroIntervention. 2025. [PMID 40887985](https://pubmed.ncbi.nlm.nih.gov/40887985/)
32. Beck S, et al. **Invasive Diagnosis of Coronary Functional Disorders Causing Angina Pectoris.** Eur Cardiol. 2021. [PMID 34276812](https://pubmed.ncbi.nlm.nih.gov/34276812/)
33. Chalikias G, et al. **Slow Coronary Flow: Pathophysiology, Clinical Implications, and Therapeutic Management.** Angiology. 2021. [PMID 33779300](https://pubmed.ncbi.nlm.nih.gov/33779300/)
34. Lee SH, et al. **Clinical Relevance of Ischemia with Nonobstructive Coronary Arteries According to Coronary Microvascular Dysfunction.** J Am Heart Assoc. 2022. [PMID 35475358](https://pubmed.ncbi.nlm.nih.gov/35475358/)

---

## 3. Indices — CFR · IMR · MRR · absolute flow

The model **distinguishes a ratio from a resistance**. The index whose interpretation survives
treatment is the hyperaemic resistance, and the references in this section supply the definitions and the limits of
those indices.

35. Fearon WF, et al. **Novel index for invasively assessing the coronary microcirculation.** Circulation. 2003. [PMID 12821539](https://pubmed.ncbi.nlm.nih.gov/12821539/) — the original definition of IMR (Pd × Tmn). The model's `IMR = MR_hyp × 8.35` scaling is set so that a normal territory reads 18 U.
36. Martínez GJ, et al. **The index of microcirculatory resistance in the physiologic assessment of the coronary microcirculation.** Coron Artery Dis. 2015. [PMID 26247265](https://pubmed.ncbi.nlm.nih.gov/26247265/)
37. Yong AS, et al. **Calculation of the index of microcirculatory resistance without coronary wedge pressure measurement in the presence of epicardial stenosis.** JACC Cardiovasc Interv. 2013. [PMID 23347861](https://pubmed.ncbi.nlm.nih.gov/23347861/)
38. De Bruyne B, et al. **Coronary thermodilution to assess flow reserve: experimental validation.** Circulation. 2001. [PMID 11673336](https://pubmed.ncbi.nlm.nih.gov/11673336/)
39. Collet C, et al. **A Systematic Approach to the Evaluation of the Coronary Microcirculation Using Bolus Thermodilution: CATH CMD.** J Soc Cardiovasc Angiogr Interv. 2024. [PMID 39131992](https://pubmed.ncbi.nlm.nih.gov/39131992/)
40. Gutiérrez-Barrios A, et al. **Continuous Thermodilution Method to Assess Coronary Flow Reserve.** Am J Cardiol. 2021. [PMID 33220317](https://pubmed.ncbi.nlm.nih.gov/33220317/)
41. Candreva A, et al. **Automation of intracoronary continuous thermodilution for absolute coronary flow and microvascular resistance measurements.** Catheter Cardiovasc Interv. 2022. [PMID 35723684](https://pubmed.ncbi.nlm.nih.gov/35723684/)
42. Belmonte M, et al. **Measuring Absolute Coronary Flow and Microvascular Resistance by Thermodilution: JACC Review Topic of the Week.** J Am Coll Cardiol. 2024. [PMID 38325996](https://pubmed.ncbi.nlm.nih.gov/38325996/)
43. Pijls NHJ, et al. **Absolute coronary blood flow measurement and the principle of microvascular resistance reserve.** Cardiovasc Interv Ther. 2026. [PMID 41432885](https://pubmed.ncbi.nlm.nih.gov/41432885/) — the principle of MRR. The model predicts that MRR **fails to correct for a drug's fall in blood pressure** (because the same pressure drops both at rest and at hyperaemia, the correction terms cancel and MRR becomes identically CFR/FFR). What MRR is of value for is epicardial stenosis, and those two problems are often confused.
44. Sakai K, et al. **Impact of vessel volume on thermodilution measurements in patients with coronary microvascular dysfunction.** Catheter Cardiovasc Interv. 2024. [PMID 38566527](https://pubmed.ncbi.nlm.nih.gov/38566527/)
45. Gallinoro E, et al. **Microvascular Dysfunction in Patients With Type II Diabetes Mellitus: Invasive Assessment of Absolute Coronary Blood Flow and Microvascular Resistance Reserve.** Front Cardiovasc Med. 2021. [PMID 34738020](https://pubmed.ncbi.nlm.nih.gov/34738020/)
46. Paolisso P, et al. **Absolute coronary flow and microvascular resistance reserve in patients with severe aortic stenosis.** Heart. 2022. [PMID 35977812](https://pubmed.ncbi.nlm.nih.gov/35977812/)
47. Johnson NP, et al. **Invasive FFR and Noninvasive CFR in the Evaluation of Ischemia: What Is the Future?** J Am Coll Cardiol. 2016. [PMID 27282899](https://pubmed.ncbi.nlm.nih.gov/27282899/)
48. Simova I. **Coronary Flow Velocity Reserve Assessment with Transthoracic Doppler Echocardiography.** Eur Cardiol. 2015. [PMID 30310417](https://pubmed.ncbi.nlm.nih.gov/30310417/)
49. Olsen RH, et al. **Coronary flow velocity reserve by echocardiography: feasibility, reproducibility and agreement with PET.** Cardiovasc Ultrasound. 2016. [PMID 27267255](https://pubmed.ncbi.nlm.nih.gov/27267255/)
50. Loftspring E, et al. **Angiography-Derived Versus Coronary Guidewire-Derived Index of Microcirculatory Resistance in Patients With INOCA.** J Soc Cardiovasc Angiogr Interv. 2025. [PMID 41324045](https://pubmed.ncbi.nlm.nih.gov/41324045/)
51. Gao B, et al. **Quantitative Flow Ratio-Derived Index of Microcirculatory Resistance as a Novel Tool to Identify Microcirculatory Function in Patients with INOCA.** Cardiology. 2024. [PMID 37839404](https://pubmed.ncbi.nlm.nih.gov/37839404/)
52. Wang S, et al. **Myocardial Blood Flow Quantification Using Stress Cardiac Magnetic Resonance Improves Detection of Coronary Artery Disease.** JACC Cardiovasc Imaging. 2024. [PMID 39297850](https://pubmed.ncbi.nlm.nih.gov/39297850/)
53. Rasmussen LD, et al. **Myocardial Blood Flow by Magnetic Resonance in Patients With Suspected Coronary Stenosis: Comparison to PET and Invasive Physiology.** Circ Cardiovasc Imaging. 2024. [PMID 38889213](https://pubmed.ncbi.nlm.nih.gov/38889213/)
54. Rasmussen LD, et al. **Impact of Absolute Myocardial Blood Flow Quantification on the Diagnostic Performance of PET-Based Perfusion Scans Using 82Rubidium.** Circ Cardiovasc Imaging. 2024. [PMID 38227687](https://pubmed.ncbi.nlm.nih.gov/38227687/)
55. Soman P, et al. **Absolute Myocardial Blood Flow Quantification With PET: Should Diagnostic Cutoffs Be Tracer Specific?** Circ Cardiovasc Imaging. 2025. [PMID 40438936](https://pubmed.ncbi.nlm.nih.gov/40438936/)
56. Klein R, et al. **Selection of PET Camera and Implications on the Reliability and Accuracy of Absolute Myocardial Blood Flow Quantification.** Curr Cardiol Rep. 2020. [PMID 32770426](https://pubmed.ncbi.nlm.nih.gov/32770426/)
57. Hagemann CE, et al. **Quantitative myocardial blood flow with Rubidium-82 PET: a clinical perspective.** Am J Nucl Med Mol Imaging. 2015. [PMID 26550537](https://pubmed.ncbi.nlm.nih.gov/26550537/)
58. Bhave NM, et al. **Considerations when measuring myocardial perfusion reserve by cardiovascular magnetic resonance using regadenoson.** J Cardiovasc Magn Reson. 2012. [PMID 23272658](https://pubmed.ncbi.nlm.nih.gov/23272658/)
59. Wöhrle J, et al. **Myocardial perfusion reserve in cardiovascular magnetic resonance: Correlation to coronary microvascular dysfunction.** J Cardiovasc Magn Reson. 2006. [PMID 17060099](https://pubmed.ncbi.nlm.nih.gov/17060099/)
60. Taqueti VR, et al. **Clinical significance of noninvasive coronary flow reserve assessment in patients with ischemic heart disease.** Curr Opin Cardiol. 2016. [PMID 27652814](https://pubmed.ncbi.nlm.nih.gov/27652814/)
61. Taqueti VR, et al. **The role of positron emission tomography in the evaluation of myocardial ischemia in women.** J Nucl Cardiol. 2016. [PMID 27488383](https://pubmed.ncbi.nlm.nih.gov/27488383/)
62. Gunasekaran V, et al. **Assessment of coronary microvascular dysfunction in INOCA using 13N-ammonia PET: Lack of correlation with angiographic flow grades.** J Nucl Cardiol. 2026. [PMID 41407151](https://pubmed.ncbi.nlm.nih.gov/41407151/)
63. Singh H, et al. **Potential Role of 13N-NH3 Cardiac PET in Monitoring Treatment Response in Patients with Microvascular Angina.** Indian J Nucl Med. 2025. [PMID 40735747](https://pubmed.ncbi.nlm.nih.gov/40735747/)
64. Zaman MU, et al. **Cardiac Positron Emission Tomography Myocardial Perfusion Imaging: Seeing Beyond Perfusion.** World J Nucl Med. 2026. [PMID 42395164](https://pubmed.ncbi.nlm.nih.gov/42395164/)

---

## 4. Prognosis

65. **Kelshiker MA, et al. Coronary flow reserve and cardiovascular outcomes: a systematic review and meta-analysis. Eur Heart J. 2022.** [PMID 34849697](https://pubmed.ncbi.nlm.nih.gov/34849697/)
    **Figures used in the model:** with an abnormal CFR, MACE HR 3.42 (95% CI 2.92–3.99). **Per 0.1 unit fall in CFR, HR for death 1.16 (1.04–1.29) and for MACE 1.08.** In patients with non-obstructive coronary arteries, the HR for death with an abnormal CFR is 5.44 (3.78–7.83). The model's hazard function `dCHMORT/dt = H0_MORT·exp(ln1.16·(2.5−CFR)/0.1)` uses this slope directly.
66. Gdowski MA, et al. **Association of Isolated Coronary Microvascular Dysfunction With Mortality and Major Adverse Cardiac Events: A Systematic Review and Meta-Analysis.** J Am Heart Assoc. 2020. [PMID 32345133](https://pubmed.ncbi.nlm.nih.gov/32345133/)
67. Luo X, et al. **Impact of Isolated Coronary Microvascular Disease Diagnosed Using Various Measurement Modalities on Prognosis: An Updated Systematic Review and Meta-Analysis.** Cardiology. 2024. [PMID 37708863](https://pubmed.ncbi.nlm.nih.gov/37708863/)
68. Gallinoro E, et al. **Prognostic Value of Microvascular Resistance Reserve in Coronary Artery Disease: A Systematic Review and Meta-Analysis.** JACC Cardiovasc Interv. 2026. [PMID 41881651](https://pubmed.ncbi.nlm.nih.gov/41881651/)
69. Seitz A, et al. **Prognostic implications of coronary artery stenosis and coronary spasm in patients with stable angina: 5-year follow-up of the ACOVA study.** Coron Artery Dis. 2020. [PMID 32168049](https://pubmed.ncbi.nlm.nih.gov/32168049/)
70. Weber BN, et al. **Impaired Coronary Vasodilator Reserve and Adverse Prognosis in Patients With Systemic Inflammatory Disorders.** JACC Cardiovasc Imaging. 2021. [PMID 33744132](https://pubmed.ncbi.nlm.nih.gov/33744132/)
71. Shah NR, et al. **Prognostic Value of Coronary Flow Reserve in Patients with Dialysis-Dependent ESRD.** J Am Soc Nephrol. 2016. [PMID 26459635](https://pubmed.ncbi.nlm.nih.gov/26459635/)
72. Huck DM, et al. **Comparative effectiveness of PET and SPECT myocardial perfusion imaging for predicting risk in patients with cardiometabolic disease.** J Nucl Cardiol. 2024. [PMID 38996910](https://pubmed.ncbi.nlm.nih.gov/38996910/)
73. Gould KL, et al. **Subendocardial and Transmural Myocardial Ischemia: Clinical Characteristics, Prevalence, and Outcomes With and Without Revascularization.** JACC Cardiovasc Imaging. 2023. [PMID 36599572](https://pubmed.ncbi.nlm.nih.gov/36599572/)
74. Eftekhari A, et al. **Changes in microvascular resistance following percutaneous coronary intervention — From the ILIAS global registry.** Int J Cardiol. 2023. [PMID 37633364](https://pubmed.ncbi.nlm.nih.gov/37633364/)
75. Jeyaprakash P, et al. **Index of Microcirculatory Resistance to predict microvascular obstruction in STEMI: systematic review and meta-analysis.** Catheter Cardiovasc Interv. 2024. [PMID 38179600](https://pubmed.ncbi.nlm.nih.gov/38179600/)
76. Zhang Y, et al. **Prognostic Value of Coronary Angiography-Derived Index of Microcirculatory Resistance in NSTEMI Patients.** JACC Cardiovasc Interv. 2024. [PMID 39115479](https://pubmed.ncbi.nlm.nih.gov/39115479/)
77. Zheng Y, et al. **Prognostic Value of Angiography-Derived IMR in Patients With Intermediate Coronary Stenosis.** JACC Cardiovasc Interv. 2025. [PMID 39880572](https://pubmed.ncbi.nlm.nih.gov/39880572/)
78. Chen D, et al. **Combined risk estimates of diabetes and angiography-derived IMR in NSTEMI.** Cardiovasc Diabetol. 2024. [PMID 39152477](https://pubmed.ncbi.nlm.nih.gov/39152477/)
79. Zhang Y, et al. **Prognostic Value of Angiography-Derived Index of Microvascular Resistance in Hypertrophic Cardiomyopathy.** MedComm. 2025. [PMID 40717902](https://pubmed.ncbi.nlm.nih.gov/40717902/)

---

## 5. Autoregulation and transmural perfusion — the physics of the model

**"Supply lives only in diastole"** is the second axis of the model, and the following references supply the DPTI/SPTI
formalism and the double entry of heart rate (into demand and into perfusion time).

80. **Duncker DJ, Bache RJ. Regulation of coronary blood flow during exercise. Physiol Rev. 2008;88:1009–86.** [PMID 18626066](https://pubmed.ncbi.nlm.nih.gov/18626066/) — the basis for the metabolic autoregulation structure of the model, a resting oxygen extraction ratio of ~70%, the ceiling on maximal extraction, and K_ATP- and adenosine-mediated dilatation.
81. Duncker DJ, et al. **Regulation of coronary resistance vessel tone in response to exercise.** J Mol Cell Cardiol. 2012. [PMID 22037538](https://pubmed.ncbi.nlm.nih.gov/22037538/)
82. Duncker DJ, et al. **Role of K+ATP channels in coronary vasodilation during exercise.** Circulation. 1993. [PMID 8353886](https://pubmed.ncbi.nlm.nih.gov/8353886/)
83. Duncker DJ, et al. **Role of K+ATP channels and adenosine in the regulation of coronary blood flow during exercise with normal and restricted coronary blood flow.** J Clin Invest. 1996. [PMID 8613554](https://pubmed.ncbi.nlm.nih.gov/8613554/)
84. Hoffman JI, Buckberg GD. **Transmural myocardial perfusion.** Prog Cardiovasc Dis. 1987. [PMID 2953043](https://pubmed.ncbi.nlm.nih.gov/2953043/) — the structure in which the subendocardium is perfused in diastole only and takes extravascular compression as a series resistance (the model's `W_ENDO`, `RCOMP_K`, `PHI_SYS`).
85. **Buckberg GD, et al. Ischemia in aortic stenosis: hemodynamic prediction. Am J Cardiol. 1975.** [PMID 1130286](https://pubmed.ncbi.nlm.nih.gov/1130286/) — the origin of the DPTI/SPTI (=SEVR) formalism.
86. Brazier JR, et al. **Effects of tachycardia on the adequacy of subendocardial oxygen delivery in experimental aortic stenosis.** Am Heart J. 1975. [PMID 1155327](https://pubmed.ncbi.nlm.nih.gov/1155327/)
87. Brazier JR, et al. **Papillary muscle ischemia with patent coronary arteries.** Surgery. 1975. [PMID 1166409](https://pubmed.ncbi.nlm.nih.gov/1166409/)
88. Canty JM Jr, et al. **Effect of tachycardia on regional function and transmural myocardial perfusion during graded coronary pressure reduction in conscious dogs.** Circulation. 1990. [PMID 2225378](https://pubmed.ncbi.nlm.nih.gov/2225378/)
89. Buck JD, et al. **Changes in ischemic blood flow distribution and dynamic severity of a coronary stenosis induced by beta blockade in the canine heart.** Circulation. 1981. [PMID 6115724](https://pubmed.ncbi.nlm.nih.gov/6115724/)
90. Buck JD, et al. **Effects of sotalol and vagal stimulation on ischemic myocardial blood flow distribution in the canine heart.** J Pharmacol Exp Ther. 1981. [PMID 7463353](https://pubmed.ncbi.nlm.nih.gov/7463353/)
91. Chemla D, et al. **Subendocardial viability ratio estimated by arterial tonometry: a critical evaluation in elderly hypertensive patients with increased aortic stiffness.** Clin Exp Pharmacol Physiol. 2008. [PMID 18346166](https://pubmed.ncbi.nlm.nih.gov/18346166/)
92. Reitan JA, et al. **A computer evaluation of the ratio of the diastolic pressure-time index to the time-tension index from three arterial sites in dogs.** J Clin Monit. 1986. [PMID 3711953](https://pubmed.ncbi.nlm.nih.gov/3711953/)
93. Kissling G, et al. **Mechanical determinants of myocardial oxygen consumption with special reference to external work and efficiency.** Cardiovasc Res. 1992. [PMID 1451165](https://pubmed.ncbi.nlm.nih.gov/1451165/) — the tension-time term of the model's MVO2 equation.
94. Balady GJ, et al. **Comparison of determinants of myocardial oxygen consumption during arm and leg exercise in normal persons.** Am J Cardiol. 1986. [PMID 3717042](https://pubmed.ncbi.nlm.nih.gov/3717042/)
95. Richalet JP, et al. **Myocardial oxygen extraction and oxygen-hemoglobin equilibrium curve during moderate exercise.** Eur J Appl Physiol. 1981. [PMID 7197622](https://pubmed.ncbi.nlm.nih.gov/7197622/) — the model's `E_REST = 0.70` and `E_MAX = 0.80`.

---

## 6. Endothelial biology — NO · BH4 · ROS · ADMA

96. Yuyun MF, et al. **Endothelial dysfunction, endothelial nitric oxide bioavailability, tetrahydrobiopterin, and 5-methyltetrahydrofolate in cardiovascular disease.** Microvasc Res. 2018. [PMID 29596860](https://pubmed.ncbi.nlm.nih.gov/29596860/)
97. Bendall JK, et al. **Tetrahydrobiopterin in cardiovascular health and disease.** Antioxid Redox Signal. 2014. [PMID 24294830](https://pubmed.ncbi.nlm.nih.gov/24294830/) — the model's BH4↔BH2 oxidation and eNOS uncoupling (`KBH4_OX`).
98. Cherng TW, et al. **Mechanisms of diesel-induced endothelial nitric oxide synthase dysfunction in coronary arterioles.** Environ Health Perspect. 2011. [PMID 20870565](https://pubmed.ncbi.nlm.nih.gov/20870565/)
99. Landim MB, et al. **Asymmetric dimethylarginine (ADMA) and endothelial dysfunction: implications for atherogenesis.** Clinics. 2009. [PMID 19488614](https://pubmed.ncbi.nlm.nih.gov/19488614/)
100. Böger RH. **Association of asymmetric dimethylarginine and endothelial dysfunction.** Clin Chem Lab Med. 2003. [PMID 14656027](https://pubmed.ncbi.nlm.nih.gov/14656027/)
101. Bełtowski J, et al. **Asymmetric dimethylarginine (ADMA) as a target for pharmacotherapy.** Pharmacol Rep. 2006. [PMID 16702618](https://pubmed.ncbi.nlm.nih.gov/16702618/)
102. Chan NN, et al. **ADMA: a potential link between endothelial dysfunction and cardiovascular diseases in insulin resistance syndrome?** Diabetologia. 2002. [PMID 12488950](https://pubmed.ncbi.nlm.nih.gov/12488950/)
103. Thengchaisri N, et al. **H2O2 Mediates VEGF- and Flow-Induced Dilations of Coronary Arterioles in Early Type 1 Diabetes: Role of Vascular Arginase and PI3K-Linked eNOS Uncoupling.** Int J Mol Sci. 2022. [PMID 36613929](https://pubmed.ncbi.nlm.nih.gov/36613929/)
104. Mahmoud AM, et al. **Nox2 contributes to hyperinsulinemia-induced redox imbalance and impaired vascular function.** Redox Biol. 2017. [PMID 28600985](https://pubmed.ncbi.nlm.nih.gov/28600985/)
105. Younis W, et al. **Soluble guanylyl cyclase, the NO receptor, drives vasorelaxation via endothelial S-nitrosation.** Proc Natl Acad Sci USA. 2025. [PMID 41037641](https://pubmed.ncbi.nlm.nih.gov/41037641/)
106. Friebe A, et al. **NO-GC in cells 'off the beaten track'.** Nitric Oxide. 2018. [PMID 29626542](https://pubmed.ncbi.nlm.nih.gov/29626542/)
107. Xiao S, et al. **Soluble Guanylate Cyclase Stimulators and Activators: Where are We and Where to Go?** Mini Rev Med Chem. 2019. [PMID 31362687](https://pubmed.ncbi.nlm.nih.gov/31362687/)
108. Torfgård KE, Ahlner J. **Mechanisms of action of nitrates.** Cardiovasc Drugs Ther. 1994. [PMID 7873467](https://pubmed.ncbi.nlm.nih.gov/7873467/)
109. Al-Badri A, et al. **Peripheral Microvascular Function Reflects Coronary Vascular Function.** Arterioscler Thromb Vasc Biol. 2019. [PMID 31018659](https://pubmed.ncbi.nlm.nih.gov/31018659/)
110. McChord J, et al. **Coronary Endothelial Dysfunction: Diagnostic Necessity or Futile Effort in Patients With Non-Obstructive Angina?** Catheter Cardiovasc Interv. 2025. [PMID 40745893](https://pubmed.ncbi.nlm.nih.gov/40745893/)

---

## 7. The endothelin axis

111. **Ford TJ, et al. Genetic dysregulation of endothelin-1 is implicated in coronary microvascular dysfunction. Eur Heart J. 2020.** [PMID 31972008](https://pubmed.ncbi.nlm.nih.gov/31972008/) — in 391 patients with angina, rs9349379-G was associated with a raised circulating ET-1 and with CMD. The basis for the model parameters `GENO` and `A_GENE_E = 0.30`, and the precision-medicine hypothesis of PRIZE.
112. Feng J, et al. **Endothelin-1-induced contractile responses of human coronary arterioles via endothelin-A receptors and PKC-alpha signaling pathways.** Surgery. 2010. [PMID 20079914](https://pubmed.ncbi.nlm.nih.gov/20079914/) — ETA-mediated arteriolar constriction (the model's `F_ETA = 0.75`, `F_ETB2 = 0.25`).
113. Dashwood MR, et al. **Regional variations in endothelin-1 and its receptor subtypes in human coronary vasculature.** Endothelium. 1998. [PMID 9832333](https://pubmed.ncbi.nlm.nih.gov/9832333/)
114. DeFily DV, et al. **Endothelin antagonists block alpha1-adrenergic constriction of coronary arterioles.** Am J Physiol. 1999. [PMID 10070088](https://pubmed.ncbi.nlm.nih.gov/10070088/) — the coupling of the α1 and ET axes (the model's `K_A1_TN`).
115. Lamping KG, et al. **Effects of 17 beta-estradiol on coronary microvascular responses to endothelin-1.** Am J Physiol. 1996. [PMID 8853349](https://pubmed.ncbi.nlm.nih.gov/8853349/)
116. Sauvageau S, et al. **Evaluation of endothelin-1-induced pulmonary vasoconstriction following myocardial infarction.** Exp Biol Med. 2006. [PMID 16741009](https://pubmed.ncbi.nlm.nih.gov/16741009/)

---

## 8. Rho-kinase and spasm

117. Shimokawa H. **2014 Williams Harvey Lecture: importance of coronary vasomotion abnormalities—from bench to bedside.** Eur Heart J. 2014. [PMID 25354517](https://pubmed.ncbi.nlm.nih.gov/25354517/)
118. Shimokawa H. **Cellular and molecular mechanisms of coronary artery spasm.** Jpn Circ J. 2000. [PMID 10651199](https://pubmed.ncbi.nlm.nih.gov/10651199/) — the model's ROCK → inhibition of MLC dephosphorylation → Ca-independent sensitisation.
119. Yoo SY, Kim JY. **Recent insights into the mechanisms of vasospastic angina.** Korean Circ J. 2009. [PMID 20049135](https://pubmed.ncbi.nlm.nih.gov/20049135/)
120. Oi K, et al. **Remnant lipoproteins from patients with sudden cardiac death enhance coronary vasospastic activity through upregulation of Rho-kinase.** Arterioscler Thromb Vasc Biol. 2004. [PMID 15044207](https://pubmed.ncbi.nlm.nih.gov/15044207/) — the pathway by which LDL and remnant lipoproteins raise `ROCK_D` in the model.
121. **Masumoto A, et al. Suppression of coronary artery spasm by the Rho-kinase inhibitor fasudil in patients with vasospastic angina. Circulation. 2002.** [PMID 11927519](https://pubmed.ncbi.nlm.nih.gov/11927519/) — the basis for the fasudil arm of the model (`EM_FAS_RK = 0.72`).
122. Otsuka T, et al. **Administration of the Rho-kinase inhibitor fasudil following nitroglycerin additionally dilates the site of coronary spasm.** Coron Artery Dis. 2008. [PMID 18300747](https://pubmed.ncbi.nlm.nih.gov/18300747/) — the observation that a component of the spasm is left unrelieved by nitrate (why nitrate acts on the epicardium only in the model).
123. Mohri M, et al. **Angina pectoris caused by coronary microvascular spasm.** Lancet. 1998. [PMID 9643687](https://pubmed.ncbi.nlm.nih.gov/9643687/)
124. Sun H, et al. **Coronary microvascular spasm causes myocardial ischemia in patients with vasospastic angina.** J Am Coll Cardiol. 2002. [PMID 11869851](https://pubmed.ncbi.nlm.nih.gov/11869851/) — microvascular spasm produces ischaemia without a stenosis → why the model writes spasm not as "tone" but as **vessel occlusion + loss of driving pressure** (bug B20).
125. Ong P, et al. **High prevalence of a pathological response to acetylcholine testing in patients with stable angina pectoris and unobstructed coronary arteries. The ACOVA Study.** J Am Coll Cardiol. 2012. [PMID 22322081](https://pubmed.ncbi.nlm.nih.gov/22322081/)
126. Seitz A, et al. **Characterization and implications of intracoronary hemodynamic assessment during coronary spasm provocation testing.** Clin Res Cardiol. 2023. [PMID 37195455](https://pubmed.ncbi.nlm.nih.gov/37195455/)
127. Feenstra RGT, et al. **Haemodynamic characterisation of different endotypes in coronary artery vasospasm in reaction to acetylcholine.** Int J Cardiol Heart Vasc. 2022. [PMID 36017267](https://pubmed.ncbi.nlm.nih.gov/36017267/)
128. Feenstra RGT, et al. **Post-spastic flow recovery time to document vasospasm induced ischemia during acetylcholine provocation testing.** Int J Cardiol Heart Vasc. 2023. [PMID 37275626](https://pubmed.ncbi.nlm.nih.gov/37275626/)
129. Feenstra RGT, et al. **Do ECG changes induced during intracoronary vasospasm provocation testing reflect those during spontaneous angina episodes in vasospastic angina?** Eur Heart J Case Rep. 2024. [PMID 39161720](https://pubmed.ncbi.nlm.nih.gov/39161720/)
130. Huang J, et al. **Invasive Evaluation for Coronary Vasospasm.** US Cardiol. 2023. [PMID 39493950](https://pubmed.ncbi.nlm.nih.gov/39493950/)
131. Aswathappa S, et al. **A Comprehensive Literature Review Discussing Diagnostic Challenges of Prinzmetal or Vasospastic Angina.** Cureus. 2025. [PMID 40486459](https://pubmed.ncbi.nlm.nih.gov/40486459/)
132. Mehta HH, et al. **The Spontaneous Coronary Slow-Flow Phenomenon: Reversal by Intracoronary Nicardipine.** J Invasive Cardiol. 2019. [PMID 30555052](https://pubmed.ncbi.nlm.nih.gov/30555052/)
133. Sykes R, et al. **Myocardial Bridging Independently Associates With Coronary Artery Spasm.** JACC Cardiovasc Interv. 2026. [PMID 42508854](https://pubmed.ncbi.nlm.nih.gov/42508854/)
134. Toya T, et al. **Coronary Endothelial Dysfunction and Vasomotor Dysregulation in Myocardial Bridging.** J Cardiovasc Dev Dis. 2025. [PMID 39997488](https://pubmed.ncbi.nlm.nih.gov/39997488/)

---

## 9. Adenosine, K_ATP, nociception

This section underpins the least obvious prediction of the model. **Because the functional endotype has a normal
minimal resistance, its subendocardial supply-demand deficit is effectively zero at any workload, and so its
angina is carried not by ischaemia but by afferent (A1-adenosine) signalling and central sensitisation.**

135. **Elliott PM, et al. Effect of oral aminophylline in patients with angina and normal coronary arteriograms (cardiac syndrome X). Heart. 1997;77:523.** [PMID 9227295](https://pubmed.ncbi.nlm.nih.gov/9227295/) — the observation the model reproduces: adenosine receptor blockade barely moves CFR yet lengthens exercise time. The model predicts that this effect is **confined to the functional endotype** and that in the structural endotype it may be harmful, because A2A blockade costs the dilator reserve.
136. Zhou X, et al. **A1 adenosine receptor negatively modulates coronary reactive hyperemia via counteracting A2A-mediated H2O2 production and KATP opening.** Am J Physiol Heart Circ Physiol. 2013. [PMID 24043252](https://pubmed.ncbi.nlm.nih.gov/24043252/) — the opposing actions of A1 and A2A (the model's `EM_AMI_A1` against `EM_AMI_A2`).
137. Peart JN, Headrick JP. **Adenosinergic cardioprotection: multiple receptors, multiple pathways.** Pharmacol Ther. 2007. [PMID 17408751](https://pubmed.ncbi.nlm.nih.gov/17408751/)
138. Riou LM, et al. **Influence of propranolol, enalaprilat, verapamil, and caffeine on adenosine A2A-receptor-mediated coronary vasodilation.** J Am Coll Cardiol. 2002. [PMID 12427424](https://pubmed.ncbi.nlm.nih.gov/12427424/) — caffeine and the methylxanthines interfere with adenosine dilatation → a conflict between the diagnostic test and the treatment.
139. Niiya K, et al. **Glibenclamide reduces the coronary vasoactivity of adenosine receptor agonists.** J Pharmacol Exp Ther. 1994. [PMID 7965706](https://pubmed.ncbi.nlm.nih.gov/7965706/)
140. Lanza GA, et al. **Effect of spinal cord stimulation on spontaneous and stress-induced angina and 'ischemia-like' ST-segment depression in patients with cardiac syndrome X.** Eur Heart J. 2005. [PMID 15642701](https://pubmed.ncbi.nlm.nih.gov/15642701/)
141. Lanza GA, et al. **Spinal cord stimulation in patients with refractory anginal pain and normal coronary arteries.** Ital Heart J. 2001. [PMID 11214698](https://pubmed.ncbi.nlm.nih.gov/11214698/)
142. Eliasson T, et al. **Spinal cord stimulation in angina pectoris with normal coronary arteriograms.** Coron Artery Dis. 1993. [PMID 8287216](https://pubmed.ncbi.nlm.nih.gov/8287216/)
143. Lanza GA, et al. **Management of microvascular angina pectoris.** Am J Cardiovasc Drugs. 2014. [PMID 24174173](https://pubmed.ncbi.nlm.nih.gov/24174173/)
144. de Silva R, et al. **Refractory angina: mechanisms and stratified treatment in obstructive and non-obstructive chronic myocardial ischaemic syndromes.** Eur Heart J. 2025. [PMID 40590516](https://pubmed.ncbi.nlm.nih.gov/40590516/)
145. Tyrer P, et al. **Cognitive behaviour therapy for non-cardiac pain in the chest (COPIC): a multicentre randomized controlled trial with economic evaluation.** BMC Psychol. 2015. [PMID 26596540](https://pubmed.ncbi.nlm.nih.gov/26596540/) — the model's `CBT` switch.
146. Eriksson-Liebon M, et al. **Long-term effects and predictors of change of internet-delivered CBT on cardiac anxiety in patients with non-cardiac chest pain: RCT.** BMC Psychiatry. 2024. [PMID 38504157](https://pubmed.ncbi.nlm.nih.gov/38504157/)
147. Thesen T, et al. **Patients with depression symptoms are more likely to experience improvements of internet-based CBT: secondary analysis in non-cardiac chest pain.** BMC Psychiatry. 2023. [PMID 37838653](https://pubmed.ncbi.nlm.nih.gov/37838653/)
148. Achem SR. **Recent developments in chest pain of undetermined origin.** Curr Gastroenterol Rep. 2000. [PMID 10957931](https://pubmed.ncbi.nlm.nih.gov/10957931/)
149. Shrestha S, Pasricha PJ. **Update on noncardiac chest pain.** Dig Dis. 2000. [PMID 11279332](https://pubmed.ncbi.nlm.nih.gov/11279332/)

---

## 10. Structural remodelling and comorbidity

150. Camici PG, et al. **Coronary microvascular dysfunction in hypertrophy and heart failure.** Cardiovasc Res. 2020. [PMID 31999329](https://pubmed.ncbi.nlm.nih.gov/31999329/) — the model's LVH and capillary density mismatch (`CAPD`, `LVH`).
151. Paulus WJ, Tschöpe C. **A novel paradigm for heart failure with preserved ejection fraction: comorbidities drive myocardial dysfunction and remodeling through coronary microvascular endothelial inflammation.** J Am Coll Cardiol. 2013. [PMID 23684677](https://pubmed.ncbi.nlm.nih.gov/23684677/)
152. **Shah SJ, et al. Prevalence and correlates of coronary microvascular dysfunction in heart failure with preserved ejection fraction: PROMIS-HFpEF. Eur Heart J. 2018.** [PMID 30165580](https://pubmed.ncbi.nlm.nih.gov/30165580/) — the coupling of CMD to filling pressure. The clinical basis for the pathway by which the rise in LVEDP on exercise in the model (`K_LVDP_W`, `K_ICF_LW`) cuts into subendocardial perfusion.
153. Sinha A, et al. **Coronary microvascular dysfunction and heart failure with preserved ejection fraction: what are the mechanistic links?** Curr Opin Cardiol. 2023. [PMID 37668191](https://pubmed.ncbi.nlm.nih.gov/37668191/)
154. Erhardsson M, et al. **Regional differences and coronary microvascular dysfunction in heart failure with preserved ejection fraction.** ESC Heart Fail. 2023. [PMID 37920127](https://pubmed.ncbi.nlm.nih.gov/37920127/)
155. Chandramouli C, et al. **Sex differences in proteomic correlates of coronary microvascular dysfunction among patients with HFpEF.** Eur J Heart Fail. 2022. [PMID 35060248](https://pubmed.ncbi.nlm.nih.gov/35060248/)
156. Venkateshvaran A, et al. **Association of epicardial adipose tissue with proteomics, coronary flow reserve, cardiac structure and function, and quality of life in HFpEF: PROMIS-HFpEF.** Eur J Heart Fail. 2022. [PMID 36196462](https://pubmed.ncbi.nlm.nih.gov/36196462/)
157. Mahmoud I, et al. **Epicardial adipose tissue differentiates in patients with and without coronary microvascular dysfunction.** Int J Obes. 2021. [PMID 34172829](https://pubmed.ncbi.nlm.nih.gov/34172829/)
158. Patel NH, et al. **Epicardial adipose tissue attenuation on CT in women with coronary microvascular dysfunction.** Atherosclerosis. 2024. [PMID 38944545](https://pubmed.ncbi.nlm.nih.gov/38944545/)
159. Agabiti-Rosei E, Rizzoni D. **[Structural and functional changes of the microcirculation in hypertension: influence of pharmacological therapy].** Drugs. 2003. [PMID 12708883](https://pubmed.ncbi.nlm.nih.gov/12708883/) — inward media/lumen remodelling and its partial reversal by an ACE inhibitor (the model's `EM_RAM_ML`, `TAU_ML`).
160. Feihl F, et al. **The macrocirculation and microcirculation of hypertension.** Curr Hypertens Rep. 2009. [PMID 19442327](https://pubmed.ncbi.nlm.nih.gov/19442327/)
161. Agabiti-Rosei E, Rizzoni D. **From macro- to microcirculation: benefits in hypertension and diabetes.** J Hypertens Suppl. 2008. [PMID 19363848](https://pubmed.ncbi.nlm.nih.gov/19363848/)
162. Sezer M, et al. **Bimodal Pattern of Coronary Microvascular Involvement in Diabetes Mellitus.** J Am Heart Assoc. 2016. [PMID 27930353](https://pubmed.ncbi.nlm.nih.gov/27930353/)
163. Niewiara Ł, et al. **Impaired coronary flow reserve in patients with poor type 2 diabetes control.** Cardiol J. 2024. [PMID 36342032](https://pubmed.ncbi.nlm.nih.gov/36342032/)
164. Huang R, et al. **Relationship between glycosylated hemoglobin A1c and coronary flow reserve in patients with type 2 diabetes.** Expert Rev Cardiovasc Ther. 2015. [PMID 25695762](https://pubmed.ncbi.nlm.nih.gov/25695762/)
165. Y-Hassan S, et al. **Coronary microvascular dysfunction in Takotsubo syndrome: cause or consequence.** Am J Cardiovasc Dis. 2021. [PMID 34084653](https://pubmed.ncbi.nlm.nih.gov/34084653/)
166. Castaldi G, et al. **Angiography-derived index of microvascular resistance in takotsubo syndrome.** Int J Cardiovasc Imaging. 2023. [PMID 36336756](https://pubmed.ncbi.nlm.nih.gov/36336756/)
167. Chitturi KR, et al. **Coronary microvascular dysfunction and cancer therapy-related cardiovascular toxicity.** Cardiovasc Revasc Med. 2024. [PMID 38789343](https://pubmed.ncbi.nlm.nih.gov/38789343/)
168. Türkoğlu C, et al. **The Relationship Between H2FPEF Score and Coronary Slow Flow Phenomenon.** Turk Kardiyol Dern Ars. 2022. [PMID 35695359](https://pubmed.ncbi.nlm.nih.gov/35695359/)
169. Amirzadegan A, et al. **Coronary slow flow phenomenon and microalbuminuria.** Turk Kardiyol Dern Ars. 2019. [PMID 31802772](https://pubmed.ncbi.nlm.nih.gov/31802772/)

---

## 11. Randomised trials and therapeutics

These are the five clinical trial anchors the model has to reproduce. Three of them are **negative**
results, and half of what the model does is to decompose, quantitatively, why those negative results
were negative.

### 11.1 Stratified therapy — CorMicA

170. **Ford TJ, et al. Stratified Medical Therapy Using Invasive Coronary Function Testing in Angina: The CorMicA Trial. J Am Coll Cardiol. 2018;72:2841–2855.** [PMID 30266608](https://pubmed.ncbi.nlm.nih.gov/30266608/)
     **Figures used in the model:** 391 patients enrolled, 206 (53.7%) with obstructive lesions on angiography and 151 (39%) non-obstructive, randomised 1:1 (76 intervention / 75 blinded control). The intervention was stratified therapy linked to guidewire-based CFR, IMR, and FFR plus ACh vasoreactivity testing. **Mean SAQ summary score at 6 months +11.7 U (95% CI 5.0–18.4, p=0.001)**, EQ-5D +0.10 (0.01–0.18), VAS +14.5 (7.8–21.3). No difference in MACE at 6 months (2.6% against 2.6%). §X of the model takes these figures as its target and produces +5.4 U (direction right, magnitude underestimated).
171. Ford TJ, et al. **Rationale and design of the BHF CorMicA stratified medicine clinical trial.** Am Heart J. 2018. [PMID 29803987](https://pubmed.ncbi.nlm.nih.gov/29803987/)
172. Ford TJ, et al. **How to Diagnose and Manage Angina Without Obstructive Coronary Artery Disease: Lessons from CorMicA.** Interv Cardiol. 2019. [PMID 31178933](https://pubmed.ncbi.nlm.nih.gov/31178933/)
173. Heggie R, et al. **Stratified medicine using invasive coronary function testing in angina: A cost-effectiveness analysis of the BHF CorMicA trial.** Int J Cardiol. 2021. [PMID 33992700](https://pubmed.ncbi.nlm.nih.gov/33992700/)

### 11.2 Ranolazine — RWISE

174. **Bairey Merz CN, et al. A randomized, placebo-controlled trial of late Na current inhibition (ranolazine) in coronary microvascular dysfunction (CMD): impact on angina and myocardial perfusion reserve. Eur Heart J. 2016;37:1504–13.** [PMID 26614823](https://pubmed.ncbi.nlm.nih.gov/26614823/)
     **Figures used in the model:** 128 patients (96% women), ranolazine 500–1000 mg twice daily for 2 weeks, double-blind crossover. **No difference in SAQ overall.** Peak heart rate on stress −3.55 bpm (p<0.001) → the model's `EM_RAN_HR = 0.052`. Correlation of the change in SAQ-7 with the change in MPRI 0.25 (p=0.005). **Only in the CFR<2.5 subgroup** were MPRI (p=0.014), angina frequency (p=0.027), and SAQ-7 (p=0.041) improved. The dilution that §VII of the model reproduces: CMD stratum +2.3 U against +1.2 U for the whole cohort (both below the MCID).
175. Hampilos KE, et al. **Myocardial biomarkers in coronary microvascular dysfunction: Response to ranolazine.** Am Heart J Plus. 2025. [PMID 40093309](https://pubmed.ncbi.nlm.nih.gov/40093309/)
176. Zhu H, et al. **Effects of the Antianginal Drugs Ranolazine, Nicorandil, and Ivabradine on Coronary Microvascular Function in Patients With Nonobstructive Coronary Artery Disease: A Meta-analysis of RCTs.** Clin Ther. 2019. [PMID 31548105](https://pubmed.ncbi.nlm.nih.gov/31548105/)
177. Patel S, et al. **Contemporary Antianginal Therapy.** Am J Cardiovasc Drugs. 2026. [PMID 40999181](https://pubmed.ncbi.nlm.nih.gov/40999181/)

### 11.3 Zibotentan — PRIZE

178. **Morrow A, et al. Zibotentan in Microvascular Angina: A Randomized, Placebo-Controlled, Crossover Trial. Circulation. 2024.** [PMID 39217504](https://pubmed.ncbi.nlm.nih.gov/39217504/)
     **Figures used in the model:** 118 patients with microvascular angina (63.5±9.2 years, 60.2% women, 21.2% with diabetes), enriched to an rs9349379-G allele frequency of 50%, zibotentan 10 mg/day for 12 weeks, sequential crossover. In the 103 with complete data, **treadmill (Bruce) duration differed by −4.26 seconds (95% CI −19.60 to +11.06, p=0.5871)**, with no improvement in any secondary index. **Zibotentan lowered blood pressure and raised circulating ET-1.** Adverse events 60.2% against 14.4% on placebo (p<0.001), predominantly fluid retention. §VI of the model produces −4.53 seconds (inside the observed confidence interval) and attributes all of it to fluid retention (−6.99 seconds). The fall in blood pressure is, if anything, worth +2.41 seconds (the reduction in demand outruns the loss of driving pressure).
179. Morrow AJ, et al. **Rationale and design of the MRC's Precision Medicine with Zibotentan in Microvascular Angina (PRIZE) trial.** Am Heart J. 2020. [PMID 32942043](https://pubmed.ncbi.nlm.nih.gov/32942043/)
180. Morrow A, et al. **Exercise treadmill testing for efficacy evaluation in randomized, controlled trials.** Am Heart J. 2026. [PMID 41687797](https://pubmed.ncbi.nlm.nih.gov/41687797/) — directly relevant to the conclusion of §VI of the model (that treadmill duration is a blunt instrument in this population).
181. Pasupathy S, et al. **Anti-Anginal Efficacy of Zibotentan in the Coronary Slow-Flow Phenomenon.** J Clin Med. 2024. [PMID 38592159](https://pubmed.ncbi.nlm.nih.gov/38592159/)

### 11.4 Intensive medical therapy — WARRIOR

182. **Pepine CJ, et al. Women's IschemiA TRial to Reduce Events In Non-ObstRuctive CAD (WARRIOR): a randomised controlled trial. Open Heart. 2026;13:e004115.** [PMID 41932694](https://pubmed.ncbi.nlm.nih.gov/41932694/)
     **Figures used in the model:** 2476 women with suspected ANOCA/INOCA, 71 centres. Intensive therapy (a high-intensity statin + ACEi/ARB + aspirin) against usual care. 421 events at 2.5 years (221 intensive / 200 usual), **primary endpoint HR 1.13 (95% CI 0.94–1.37, p=0.20)**, with no difference in the secondary endpoints either. **Admission for angina was the main component of MACE.** A sensitivity analysis correcting for contamination gave HR 0.74 (0.352–1.558, p=0.43). Enrolment fell short of plan, giving an older population (mean 64 years) with good blood pressure and LDL and high rates of statin and ACEi/ARB use. §IX of the model gives a true HR of 0.90 and shows that at 80% background contamination the observed HR becomes 0.98.
183. Lakshmanan S, et al. **Comparison of risk profiles of participants in the WARRIOR trial, using CCTA vs invasive coronary angiography.** Prog Cardiovasc Dis. 2024. [PMID 38547955](https://pubmed.ncbi.nlm.nih.gov/38547955/)

### 11.5 Heart rate control · calcium channel blockers · other drugs

184. Camici PG, et al. **Ivabradine in chronic stable angina: Effects by and beyond heart rate reduction.** Int J Cardiol. 2016. [PMID 27104917](https://pubmed.ncbi.nlm.nih.gov/27104917/) — corresponds directly to §III of the model (43% of the benefit of a lower heart rate is diastolic perfusion time).
185. Heusch G. **Ivabradine: Cardioprotection By and Beyond Heart Rate Reduction.** Drugs. 2016. [PMID 27041289](https://pubmed.ncbi.nlm.nih.gov/27041289/)
186. Giavarini A, et al. **The Role of Ivabradine in the Management of Angina Pectoris.** Cardiovasc Drugs Ther. 2016. [PMID 27475447](https://pubmed.ncbi.nlm.nih.gov/27475447/)
187. Bucchi A, et al. **Heart rate reduction via selective 'funny' channel blockers.** Curr Opin Pharmacol. 2007. [PMID 17267284](https://pubmed.ncbi.nlm.nih.gov/17267284/)
188. Borer JS, et al. **Characterization of the heart rate-lowering action of ivabradine, a selective I(f) current inhibitor.** Am J Ther. 2008. [PMID 18806523](https://pubmed.ncbi.nlm.nih.gov/18806523/) — the model's ivabradine PK/PD.
189. Chaudhary R, et al. **Ivabradine: Heart Failure and Beyond.** J Cardiovasc Pharmacol Ther. 2016. [PMID 26721645](https://pubmed.ncbi.nlm.nih.gov/26721645/)
190. Doesch AO, et al. **Heart rate reduction after heart transplantation with beta-blocker versus the selective If channel antagonist ivabradine.** Transplantation. 2007. [PMID 17989604](https://pubmed.ncbi.nlm.nih.gov/17989604/)
191. Rognoni A, et al. **Ivabradine: cardiovascular effects.** Recent Pat Cardiovasc Drug Discov. 2009. [PMID 19149708](https://pubmed.ncbi.nlm.nih.gov/19149708/)
192. Chen JW, et al. **Effects of short-term treatment of nicorandil on exercise-induced myocardial ischemia and abnormal cardiac autonomic activity in microvascular angina.** Am J Cardiol. 1997. [PMID 9205016](https://pubmed.ncbi.nlm.nih.gov/9205016/)
193. Hirohata A, et al. **Nicorandil prevents microvascular dysfunction resulting from PCI in patients with stable angina pectoris: a randomised study.** EuroIntervention. 2014. [PMID 24457276](https://pubmed.ncbi.nlm.nih.gov/24457276/)
194. Zhang Y, et al. **The effectiveness and safety of nicorandil in the treatment of patients with microvascular angina: protocol for systematic review and meta-analysis.** Medicine. 2021. [PMID 33466132](https://pubmed.ncbi.nlm.nih.gov/33466132/)
195. Jia Q, et al. **The effect of nicorandil in patients with cardiac syndrome X: a meta-analysis of RCTs.** Medicine. 2020. [PMID 32925783](https://pubmed.ncbi.nlm.nih.gov/32925783/)
196. Pavão RB, et al. **Aspirin plus verapamil relieves angina and perfusion abnormalities in patients with coronary microvascular dysfunction and Chagas disease.** Rev Soc Bras Med Trop. 2021. [PMID 34787258](https://pubmed.ncbi.nlm.nih.gov/34787258/)
197. Denardo SJ, et al. **Effect of phosphodiesterase type 5 inhibition on microvascular coronary dysfunction in women: a WISE ancillary study.** Clin Cardiol. 2011. [PMID 21780138](https://pubmed.ncbi.nlm.nih.gov/21780138/) — the sildenafil arm of the model.
198. Guarini G, et al. **Trimetazidine and Other Metabolic Modifiers.** Eur Cardiol. 2018. [PMID 30697354](https://pubmed.ncbi.nlm.nih.gov/30697354/) — the model's `TMZ` switch (the same ATP on less oxygen).
199. Nalbantgil S, et al. **The Effect of Trimetazidine in the Treatment of Microvascular Angina.** Int J Angiol. 1999. [PMID 9826407](https://pubmed.ncbi.nlm.nih.gov/9826407/)
200. Lanza GA. **[Therapy of microvascular angina].** Cardiologia. 1993. [PMID 7912650](https://pubmed.ncbi.nlm.nih.gov/7912650/)
201. Ferrara L, et al. **[Syndrome X and microvascular angina].** Minerva Cardioangiol. 1998. [PMID 9882962](https://pubmed.ncbi.nlm.nih.gov/9882962/)
202. Rusticali G, et al. **[The noninvasive identification of patients with angina and normal coronary arteries].** G Ital Cardiol. 1995. [PMID 8529853](https://pubmed.ncbi.nlm.nih.gov/8529853/)
203. Hinoi T, et al. **Acute effect of atorvastatin on coronary circulation measured by transthoracic Doppler echocardiography in patients without coronary artery disease by angiography.** Am J Cardiol. 2005. [PMID 15979441](https://pubmed.ncbi.nlm.nih.gov/15979441/)
204. Bountioukos M, et al. **Effect of atorvastatin on myocardial contractile reserve assessed by tissue Doppler imaging in moderately hypercholesterolemic patients without heart disease.** Am J Cardiol. 2003. [PMID 12943890](https://pubmed.ncbi.nlm.nih.gov/12943890/)
205. Mortensen MB, et al. **Influence of intensive lipid-lowering on CT-derived fractional flow reserve in patients with stable chest pain: FLOWPROMOTE.** Clin Cardiol. 2022. [PMID 36056636](https://pubmed.ncbi.nlm.nih.gov/36056636/)

---

## 12. Non-pharmacological and device therapy

206. Kissel CK, Nikoletou D. **Cardiac Rehabilitation and Exercise Prescription in Symptomatic Patients with Non-Obstructive Coronary Artery Disease—a Systematic Review.** Curr Treat Options Cardiovasc Med. 2018. [PMID 30121850](https://pubmed.ncbi.nlm.nih.gov/30121850/) — the model's `REHAB` switch.
207. Hausvater A, et al. **Cardiac Rehabilitation for Patients With INOCA and MINOCA: A Review.** J Cardiopulm Rehabil Prev. 2025. [PMID 40476778](https://pubmed.ncbi.nlm.nih.gov/40476778/)
208. Carvalho EE, et al. **Improved endothelial function and reversal of myocardial perfusion defects after aerobic physical training in a patient with microvascular myocardial ischemia.** Am J Phys Med Rehabil. 2011. [PMID 20531160](https://pubmed.ncbi.nlm.nih.gov/20531160/)
209. Akyuz A. **Exercise and Coronary Heart Disease.** Adv Exp Med Biol. 2020. [PMID 32342457](https://pubmed.ncbi.nlm.nih.gov/32342457/)
210. Hong J, et al. **Exercise training mitigates ER stress and UCP2 deficiency-associated coronary vascular dysfunction in atherosclerosis.** Sci Rep. 2021. [PMID 34326395](https://pubmed.ncbi.nlm.nih.gov/34326395/)
211. Tryon D, et al. **Coronary Sinus Reducer Improves Angina, Quality of Life, and Coronary Flow Reserve in Microvascular Dysfunction.** JACC Cardiovasc Interv. 2024. [PMID 39520443](https://pubmed.ncbi.nlm.nih.gov/39520443/) — the idea of raising coronary venous pressure to redistribute subendocardial perfusion. The model's `PV` and its subendocardial driving pressure term carry this mechanism.
212. Tebaldi M, et al. **Coronary Sinus Narrowing Improves Coronary Microcirculation Function in Patients With Refractory Angina: INROAD.** Circ Cardiovasc Interv. 2024. [PMID 38227697](https://pubmed.ncbi.nlm.nih.gov/38227697/)
213. Konigstein M, et al. **Coronary Sinus Narrowing for the Treatment of Patients With Angina and Evidence of Microvascular Dysfunction.** Can J Cardiol. 2026. [PMID 42309353](https://pubmed.ncbi.nlm.nih.gov/42309353/)
214. Tomaniak M, et al. **Coronary Sinus Reduction for Refractory Angina Caused by Microvascular Dysfunction—A Systematic Review.** J Clin Med. 2025. [PMID 41517541](https://pubmed.ncbi.nlm.nih.gov/41517541/)
215. Ashokprabhu ND, et al. **Enhanced External Counterpulsation for the Treatment of Angina With Nonobstructive Coronary Artery Disease.** Am J Cardiol. 2024. [PMID 37890564](https://pubmed.ncbi.nlm.nih.gov/37890564/)
216. Kronhaus KD, Lawson WE. **Enhanced external counterpulsation is an effective treatment for Syndrome X.** Int J Cardiol. 2009. [PMID 18590931](https://pubmed.ncbi.nlm.nih.gov/18590931/)
217. Bondesson SM, et al. **Reduced peripheral vascular reactivity in refractory angina pectoris: Effect of enhanced external counterpulsation.** J Geriatr Cardiol. 2011. [PMID 22783308](https://pubmed.ncbi.nlm.nih.gov/22783308/)

---

## 13. Sex · menopause · population

218. Waheed N, et al. **Sex differences in non-obstructive coronary artery disease.** Cardiovasc Res. 2020. [PMID 31958135](https://pubmed.ncbi.nlm.nih.gov/31958135/)
219. Jansen TPJ, et al. **Sex Differences in Coronary Function Test Results in Patients With Angina and Nonobstructive Disease.** Front Cardiovasc Med. 2021. [PMID 34722680](https://pubmed.ncbi.nlm.nih.gov/34722680/)
220. Steinberg RR, et al. **Coronary microvascular disease in women: epidemiology, mechanisms, evaluation, and treatment.** Can J Physiol Pharmacol. 2024. [PMID 38728748](https://pubmed.ncbi.nlm.nih.gov/38728748/)
221. Mathew D, et al. **Coronary microvascular dysfunction in menopausal women.** Heart. 2026. [PMID 42331611](https://pubmed.ncbi.nlm.nih.gov/42331611/) — the model's `RF_MENO` (loss of oestrogen → eNOS).
222. Shufelt CL, et al. **Sex-Specific Physiology and Cardiovascular Disease.** Adv Exp Med Biol. 2018. [PMID 30051400](https://pubmed.ncbi.nlm.nih.gov/30051400/)
223. SenthilKumar G, et al. **17β-Estradiol promotes sex-specific dysfunction in isolated human arterioles.** Am J Physiol Heart Circ Physiol. 2023. [PMID 36607795](https://pubmed.ncbi.nlm.nih.gov/36607795/)
224. Lam CSP, et al. **Sex differences in heart failure.** Eur Heart J. 2019. [PMID 31800034](https://pubmed.ncbi.nlm.nih.gov/31800034/)
225. Carlini NA, et al. **Vascular function in women with heart failure with preserved ejection fraction: a mismatch beyond diastole.** J Appl Physiol. 2025. [PMID 40839391](https://pubmed.ncbi.nlm.nih.gov/40839391/)
226. Ryk-Adamska M, et al. **Ophthalmological Microvascular Changes in ANOCA/INOCA Disease and Ophthalmological Methods to Detect Them—A Systematic Review.** J Clin Med. 2026. [PMID 41753032](https://pubmed.ncbi.nlm.nih.gov/41753032/)

---

## 14. Quality-of-life scales and trial endpoints

227. Spertus JA, et al. **Minimally Important Kansas City Cardiomyopathy Questionnaire Changes Across the Spectrum of Heart Failure Severity.** JACC Heart Fail. 2025. [PMID 40908080](https://pubmed.ncbi.nlm.nih.gov/40908080/) — consulted as the basis for setting the model's `MCID_SAQ = 10 U`. **Note: this paper is about the KCCQ, not the SAQ.** The MCID of the SAQ summary score could not be confirmed directly in this session, so the model takes the +11.7 U that CorMicA achieved as a practical standard of clinical importance and uses 10 U as the threshold. This value affects the direction of the model's conclusion (that the subgroup effect falls below the threshold), but the computed figure of +2.3 U itself is independent of the choice of threshold.
228. Cao H, et al. **Use of comparative effectiveness research for similar Chinese patent medicine for angina pectoris: a new approach based on patient-important outcomes.** Trials. 2014. [PMID 24641790](https://pubmed.ncbi.nlm.nih.gov/24641790/)
229. Picano E, et al. **The clinical use of stress echocardiography in chronic coronary syndromes and beyond coronary artery disease: EACVI clinical consensus statement.** Eur Heart J Cardiovasc Imaging. 2024. [PMID 37798126](https://pubmed.ncbi.nlm.nih.gov/37798126/)

---

## 15. The model's falsifiable predictions, and where to look

These are predictions the model puts forward without having been able to confirm them in the literature. For each one,
a measurement that could verify or refute it is given alongside.

| # | Prediction | How to test it | Supporting references |
|---|------|-----------|-----------|
| P1 | The **resting myocardial oxygen extraction ratio of the functional endotype is low, about 56% of normal**, and coronary venous oxygen saturation is correspondingly high | a single coronary sinus blood sample | pure deduction from §5 (80, 95). The rise in resting flow itself is measured in 20 |
| P2 | Reading a hyperaemic microvascular resistance ≥2.5 as "structural" **fails to distinguish remodelling from adenosine-resistant contractile tone**. Repeat the measurement after ROCK inhibition and some fall below the cut-off | repeat the hyperaemic measurement after acute Rho-kinase or ETA blockade | 121, 122, 124 (fasudil relieves the residual contractile component) |
| P3 | The benefit of ivabradine/β-blockade is greatest in the patients with **the highest resting heart rate, not the lowest CFR** | stratified randomisation, or a post hoc test for interaction | §III, 184–186 |
| P4 | Aminophylline is effective **only in the functional endotype** and may be harmful in the structural endotype through A2A blockade | an endotype-stratified crossover trial | 135, 136, 138 |
| P5 | The negative zibotentan result is not target failure but **fluid retention**. With a co-administered diuretic, or at a lower dose, the exercise time signal revives | add a diuretic co-treatment arm to the PRIZE design | 178 |
| P6 | The **hyperaemic microvascular resistance**, not CFR, should be prespecified as the treatment-response endpoint. A drug that lowers blood pressure lowers CFR but does not lower the hyperaemic resistance | re-analysis of existing trial data (trials that reported both indices) | 43, 178 |

---

## 16. Methods and tooling

- mrgsolve: <https://mrgsolve.org/>
- Graphviz: <https://graphviz.org/>
- PubMed E-utilities: <https://www.ncbi.nlm.nih.gov/books/NBK25501/>
- `cmd_reference_model.py` in this directory is a dependency-free Python reference implementation, written so that
  every equation could be run first without an R runtime; it records the **25
  real defects (B1–B25)** that the numerical work exposed in a BUG LOG at the top of the file, with a comment on
  the line each defect was on. `cmd_reference_output.txt` is the origin of every computed figure.

---

*Every PMID in this document was looked up directly with the PubMed E-utilities and its title, journal, and
year confirmed. For the quantitative anchors (20, 65, 170, 174, 178, 182) the full abstract was retrieved and
the figures checked. This is a model for educational and research purposes and cannot be used for clinical
decision-making.*
