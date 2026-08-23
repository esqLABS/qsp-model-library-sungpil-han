# Acute Respiratory Distress Syndrome (ARDS) — References
## Acute Respiratory Distress Syndrome — Annotated Reference List

This list is organised by section so that the provenance of the 17 clusters of `ards_qsp_model.dot` and of the 38 compartments and parameters of
`ards_mrgsolve_model.R` can be traced. Every entry includes a PubMed link.

> **Notation**: 🔑 = a reference that directly determined a particular parameter or structure of the model.

---

## 1. Definition · epidemiology · autopsy pathology

1. 🔑 ARDS Definition Task Force. **Acute respiratory distress syndrome: the Berlin Definition.** *JAMA* 2012;307:2526-33. — the mild/moderate/severe classification by PaO2/FiO2 ≤300/≤200/≤100. The model's `Berlin_class` output. <https://pubmed.ncbi.nlm.nih.gov/22797452/>
2. Matthay MA, et al. **A new global definition of acute respiratory distress syndrome.** *Am J Respir Crit Care Med* 2024;209:37-47. — the Kigali/HFNC extension of the definition. <https://pubmed.ncbi.nlm.nih.gov/37487152/>
3. 🔑 Bellani G, et al. **Epidemiology, patterns of care, and mortality for patients with ARDS in intensive care units in 50 countries (LUNG SAFE).** *JAMA* 2016;315:788-800. — hospital mortality of 40-46% in moderate-to-severe disease; the reference epidemiology for calibrating the model's hazard block. <https://pubmed.ncbi.nlm.nih.gov/26903337/>
4. Ashbaugh DG, et al. **Acute respiratory distress in adults.** *Lancet* 1967;2:319-23. — the first description. <https://pubmed.ncbi.nlm.nih.gov/4143721/>
5. Thille AW, et al. **Comparison of the Berlin definition for ARDS with autopsy: diffuse alveolar damage.** *Am J Respir Crit Care Med* 2013;187:761-7. — the mismatch between the clinical definition and the histological finding of DAD. <https://pubmed.ncbi.nlm.nih.gov/23370917/>
6. Matthay MA, Zemans RL, et al. **Acute respiratory distress syndrome.** *Nat Rev Dis Primers* 2019;5:18. — a mechanistic review; the skeleton of the whole structure of the model. <https://pubmed.ncbi.nlm.nih.gov/30872586/>
7. Meyer NJ, Gattinoni L, Calfee CS. **Acute respiratory distress syndrome.** *Lancet* 2021;398:622-37. <https://pubmed.ncbi.nlm.nih.gov/34217425/>

## 2. Alveolar epithelial injury and fluid clearance

8. 🔑 Ware LB, Matthay MA. **Alveolar fluid clearance is impaired in the majority of patients with acute lung injury.** *Am J Respir Crit Care Med* 2001;163:1376-83. — AFC <3%/h = raised mortality. The basis for the model's `AFC0` and `AFC_per_h` outputs. <https://pubmed.ncbi.nlm.nih.gov/11371404/>
9. Matthay MA, Folkesson HG, Clerici C. **Lung epithelial fluid transport and the resolution of pulmonary edema.** *Physiol Rev* 2002;82:569-600. — the ENaC/Na-K-ATPase mechanism. <https://pubmed.ncbi.nlm.nih.gov/12087129/>
10. 🔑 Ware LB, Matthay MA. **The acute respiratory distress syndrome.** *N Engl J Med* 2000;342:1334-49. — the BAL/plasma protein ratio >0.7 criterion for permeability oedema. <https://pubmed.ncbi.nlm.nih.gov/10793167/>
11. Zemans RL, et al. **Neutrophil transmigration triggers repair of the lung epithelium via β-catenin signaling.** *Proc Natl Acad Sci USA* 2011;108:15990-5. <https://pubmed.ncbi.nlm.nih.gov/21880956/>
12. Barkauskas CE, et al. **Type 2 alveolar cells are stem cells in adult lung.** *J Clin Invest* 2013;123:3025-36. — AT2 → AT1 regeneration; the model's `AT2P` compartment. <https://pubmed.ncbi.nlm.nih.gov/23921127/>
13. 🔑 Kobayashi Y, et al. **Persistence of a regeneration-associated, transitional alveolar epithelial cell state in pulmonary fibrosis.** *Nat Cell Biol* 2020;22:934-46. — arrest in a KRT5+/basaloid transitional state; the model's `KTGF_AT2` term. <https://pubmed.ncbi.nlm.nih.gov/32661339/>
14. Calfee CS, et al. **Plasma receptor for advanced glycation end-products and clinical outcomes in ARDS.** *Thorax* 2008;63:1083-9. — RAGE as a marker of epithelial injury. <https://pubmed.ncbi.nlm.nih.gov/18566109/>
15. Determann RM, et al. **Plasma levels of surfactant protein D and KL-6 for evaluation of lung injury.** *BMC Pulm Med* 2010;10:6. <https://pubmed.ncbi.nlm.nih.gov/20158912/>
16. Rokkam D, et al. **Claudin-4 levels are associated with intact alveolar fluid clearance in human lungs.** *Am J Pathol* 2011;179:1081-7. <https://pubmed.ncbi.nlm.nih.gov/21763677/>

## 3. Surfactant

17. Gregory TJ, et al. **Surfactant chemical composition and biophysical activity in ARDS.** *J Clin Invest* 1991;88:1976-81. — surfactant inactivation. The model's `SURF` and `KSI_PROT`. <https://pubmed.ncbi.nlm.nih.gov/1752956/>
18. Spragg RG, et al. **Effect of recombinant surfactant protein C-based surfactant on ARDS.** *N Engl J Med* 2004;351:884-92. — a negative result in adults. <https://pubmed.ncbi.nlm.nih.gov/15329426/>
19. Willson DF, et al. **Effect of exogenous surfactant (calfactant) in pediatric acute lung injury.** *JAMA* 2005;293:470-6. <https://pubmed.ncbi.nlm.nih.gov/15671432/>

## 4. Endothelium and glycocalyx

20. 🔑 Parikh SM, et al. **Excess circulating angiopoietin-2 may contribute to pulmonary vascular leak in sepsis.** *PLoS Med* 2006;3:e46. — the Ang-2/Tie2 axis; the model's `ANG2` compartment. <https://pubmed.ncbi.nlm.nih.gov/16417407/>
21. Calfee CS, et al. **Distinct molecular phenotypes of direct vs indirect ARDS in single and multi-center studies.** *Chest* 2015;147:1539-48. — differences in biomarkers (Ang-2, RAGE) between direct and indirect injury; the model's `DIRECT` parameter. <https://pubmed.ncbi.nlm.nih.gov/26033126/>
22. Schmidt EP, et al. **The pulmonary endothelial glycocalyx regulates neutrophil adhesion and lung injury during experimental sepsis.** *Nat Med* 2012;18:1217-23. <https://pubmed.ncbi.nlm.nih.gov/22820644/>
23. Ware LB, et al. **Significance of von Willebrand factor in septic and nonseptic patients with acute lung injury.** *Am J Respir Crit Care Med* 2004;170:766-72. <https://pubmed.ncbi.nlm.nih.gov/15201135/>
24. Mehta D, Malik AB. **Signaling mechanisms regulating endothelial permeability.** *Physiol Rev* 2006;86:279-367. — the RhoA/Rac1/S1P balance. <https://pubmed.ncbi.nlm.nih.gov/16371600/>
25. Garcia JG, et al. **Sphingosine 1-phosphate promotes endothelial cell barrier integrity by Edg-dependent cytoskeletal rearrangement.** *J Clin Invest* 2001;108:689-701. <https://pubmed.ncbi.nlm.nih.gov/11544274/>

## 5. Innate immunity, neutrophils, NETs

26. Grommes J, Soehnlein O. **Contribution of neutrophils to acute lung injury.** *Mol Med* 2011;17:293-307. <https://pubmed.ncbi.nlm.nih.gov/21046059/>
27. 🔑 Imai Y, et al. **Identification of oxidative stress and Toll-like receptor 4 signaling as a key pathway of acute lung injury.** *Cell* 2008;133:235-49. — the TLR4-oxidised phospholipid axis. <https://pubmed.ncbi.nlm.nih.gov/18423196/>
28. Dolinay T, et al. **Inflammasome-regulated cytokines are critical mediators of acute lung injury.** *Am J Respir Crit Care Med* 2012;185:1225-34. — NLRP3/IL-18. The model's `KNLRP`. <https://pubmed.ncbi.nlm.nih.gov/22461369/>
29. Grailer JJ, et al. **Critical role for the NLRP3 inflammasome during acute lung injury.** *J Immunol* 2014;192:5974-83. <https://pubmed.ncbi.nlm.nih.gov/24795455/>
30. Lefrançais E, et al. **Maladaptive role of neutrophil extracellular traps in pathogen-induced lung injury.** *JCI Insight* 2018;3:e98178. — the model's `NETS` compartment. <https://pubmed.ncbi.nlm.nih.gov/29415887/>
31. Middleton EA, et al. **Neutrophil extracellular traps contribute to immunothrombosis in COVID-19 ARDS.** *Blood* 2020;136:1169-79. <https://pubmed.ncbi.nlm.nih.gov/32597954/>
32. Abraham E, et al. **Efficacy and safety of LY315920Na/S-5920, a selective inhibitor of 14-kDa group IIA sPLA2, in severe sepsis.** *Crit Care Med* 2003;31:718-28. <https://pubmed.ncbi.nlm.nih.gov/12626977/>
33. **Sivelestat (STRIVE) — Zeiher BG, et al. Neutrophil elastase inhibition in acute lung injury: results of the STRIVE study.** *Crit Care Med* 2004;32:1695-702. — stopped for worse long-term mortality. <https://pubmed.ncbi.nlm.nih.gov/15286546/>
34. Aggarwal NR, King LS, D'Alessio FR. **Diverse macrophage populations mediate acute lung inflammation and resolution.** *Am J Physiol Lung Cell Mol Physiol* 2014;306:L709-25. — the M1/M2 programmes; the model's `M1`/`M2`. <https://pubmed.ncbi.nlm.nih.gov/24508730/>

## 6. The cytokine network

35. 🔑 Meduri GU, et al. **Persistent elevation of inflammatory cytokines predicts a poor outcome in ARDS.** *Chest* 1995;107:1062-73. — persistently raised cytokines = death. <https://pubmed.ncbi.nlm.nih.gov/7705118/>
36. Park WY, et al. **Cytokine balance in the lungs of patients with acute respiratory distress syndrome.** *Am J Respir Crit Care Med* 2001;164:1896-903. — IL-1β carries the greatest bioactivity in BAL. <https://pubmed.ncbi.nlm.nih.gov/11734443/>
37. Parsons PE, et al. **Lower tidal volume ventilation and plasma cytokine markers of inflammation in patients with acute lung injury.** *Crit Care Med* 2005;33:1-6. — LTVV lowers plasma IL-6 = clinical evidence of biotrauma. <https://pubmed.ncbi.nlm.nih.gov/15644641/>
38. Ranieri VM, et al. **Effect of mechanical ventilation on inflammatory mediators in patients with ARDS: a randomized controlled trial.** *JAMA* 1999;282:54-61. — 🔑 the demonstration of the biotrauma concept in humans. <https://pubmed.ncbi.nlm.nih.gov/10404912/>
39. Donnelly SC, et al. **Interleukin-8 and development of adult respiratory distress syndrome in at-risk patient groups.** *Lancet* 1993;341:643-7. <https://pubmed.ncbi.nlm.nih.gov/8095568/>

## 7. Alveolar coagulation and fibrin

40. 🔑 Ware LB, et al. **Pathogenetic and prognostic significance of altered coagulation and fibrinolysis in acute lung injury/ARDS.** *Crit Care Med* 2007;35:1821-8. — PAI-1 as an independent predictor of death; the model's `PAI1` and `FIBA`. <https://pubmed.ncbi.nlm.nih.gov/17667242/>
41. Idell S. **Coagulation, fibrinolysis, and fibrin deposition in acute lung injury.** *Crit Care Med* 2003;31(4 Suppl):S213-20. <https://pubmed.ncbi.nlm.nih.gov/12682443/>
42. Ackermann M, et al. **Pulmonary vascular endothelialitis, thrombosis, and angiogenesis in Covid-19.** *N Engl J Med* 2020;383:120-8. <https://pubmed.ncbi.nlm.nih.gov/32437596/>
43. **ATTACC/ACTIV-4a/REMAP-CAP Investigators. Therapeutic anticoagulation with heparin in critically ill patients with Covid-19.** *N Engl J Med* 2021;385:777-89. — no benefit in the critically ill. <https://pubmed.ncbi.nlm.nih.gov/34351722/>

## 8. Pulmonary oedema · lung water and the Starling balance

44. 🔑 Jozwiak M, et al. **Extravascular lung water is an independent prognostic factor in patients with ARDS.** *Crit Care Med* 2013;41:472-80. — EVLWI; the model's `EVLW` compartment and its hazard term. <https://pubmed.ncbi.nlm.nih.gov/23263578/>
45. Sakka SG, et al. **Prognostic value of extravascular lung water in critically ill patients.** *Chest* 2002;122:2080-6. <https://pubmed.ncbi.nlm.nih.gov/12475851/>
46. 🔑 **National Heart, Lung, and Blood Institute ARDS Clinical Trials Network (FACTT). Comparison of two fluid-management strategies in acute lung injury.** *N Engl J Med* 2006;354:2564-75. — conservative fluids: VFD +2.5 days, no difference in mortality. The validation target of scenario 7 of the model. <https://pubmed.ncbi.nlm.nih.gov/16714767/>
47. Grissom CK, et al. **Fluid management with a simplified conservative protocol for ARDS.** *Crit Care Med* 2015;43:288-95. <https://pubmed.ncbi.nlm.nih.gov/25599463/>

## 9. Respiratory mechanics and the baby lung

48. 🔑 Gattinoni L, Pesenti A. **The concept of "baby lung".** *Intensive Care Med* 2005;31:776-84. — the source of the model's `Aerated_frac` and `CRS0` structure. <https://pubmed.ncbi.nlm.nih.gov/15812622/>
49. 🔑 Amato MB, et al. **Driving pressure and survival in the acute respiratory distress syndrome.** *N Engl J Med* 2015;372:747-55. — a mediation analysis in 3,562 patients; ΔP is the strongest ventilator variable. The model's hazard `B_DP`. <https://pubmed.ncbi.nlm.nih.gov/25693014/>
50. 🔑 Gattinoni L, et al. **Ventilator-related causes of lung injury: the mechanical power.** *Intensive Care Med* 2016;42:1567-75. — the MP formula and the threshold of ~12 J/min; the model's `MechPower_Jmin` and `B_MP`. <https://pubmed.ncbi.nlm.nih.gov/27620287/>
51. Serpa Neto A, et al. **Mechanical power of ventilation is associated with mortality in critically ill patients.** *Intensive Care Med* 2018;44:1914-22. <https://pubmed.ncbi.nlm.nih.gov/30291376/>
52. Chiumello D, et al. **Lung stress and strain during mechanical ventilation for ARDS.** *Am J Respir Crit Care Med* 2008;178:346-55. — a strain threshold of ~1.5-2.0. The model's `STRAIN_TH`. <https://pubmed.ncbi.nlm.nih.gov/18451319/>
53. Cressoni M, et al. **Lung inhomogeneity in patients with acute respiratory distress syndrome.** *Am J Respir Crit Care Med* 2014;189:149-58. — stress raisers ×2; the mechanism behind the model's `PRONE_HOM`. <https://pubmed.ncbi.nlm.nih.gov/24261322/>
54. Gattinoni L, et al. **Lung recruitment in patients with the acute respiratory distress syndrome.** *N Engl J Med* 2006;354:1775-86. — the model's `RECRUIT` parameter. <https://pubmed.ncbi.nlm.nih.gov/16641394/>
55. Chen L, et al. **Potential for lung recruitment estimated by the recruitment-to-inflation ratio.** *Am J Respir Crit Care Med* 2020;201:178-87. <https://pubmed.ncbi.nlm.nih.gov/31577153/>

## 10. VILI·P-SILI (Ventilator- and Patient-Self-Inflicted Lung Injury)

56. Dreyfuss D, Saumon G. **Ventilator-induced lung injury: lessons from experimental studies.** *Am J Respir Crit Care Med* 1998;157:294-323. — volutrauma vs barotrauma. <https://pubmed.ncbi.nlm.nih.gov/9445314/>
57. Slutsky AS, Ranieri VM. **Ventilator-induced lung injury.** *N Engl J Med* 2013;369:2126-36. — a review of biotrauma; the model's `KBIO`. <https://pubmed.ncbi.nlm.nih.gov/24283226/>
58. 🔑 Brochard L, Slutsky A, Pesenti A. **Mechanical ventilation to minimize progression of lung injury in acute respiratory failure.** *Am J Respir Crit Care Med* 2017;195:438-42. — the P-SILI concept; the model's `Effort_PSILI` and `EFF_GAIN`. <https://pubmed.ncbi.nlm.nih.gov/27626833/>
59. Yoshida T, et al. **Spontaneous effort causes occult pendelluft during mechanical ventilation.** *Am J Respir Crit Care Med* 2013;188:1420-7. <https://pubmed.ncbi.nlm.nih.gov/24199628/>
60. Barnes T, et al. **Mechanical ventilation for the acute respiratory distress syndrome: hyperoxia and absorption atelectasis.** *Respir Care* 2018;63:1442-51 (a related review). <https://pubmed.ncbi.nlm.nih.gov/30154129/>
61. **ICU-ROX Investigators. Conservative oxygen therapy during mechanical ventilation in the ICU.** *N Engl J Med* 2020;382:989-98. <https://pubmed.ncbi.nlm.nih.gov/31613432/>

## 11. Gas exchange and dead space

62. 🔑 Nuckton TJ, et al. **Pulmonary dead-space fraction as a risk factor for death in the acute respiratory distress syndrome.** *N Engl J Med* 2002;346:1281-6. — Vd/Vt as an independent risk of death; the model's hazard `B_VD`. <https://pubmed.ncbi.nlm.nih.gov/11973365/>
63. Dantzker DR, et al. **Ventilation-perfusion distributions in the adult respiratory distress syndrome.** *Am Rev Respir Dis* 1979;120:1039-52. — the Riley shunt analysis; the basis for calibrating the model's `KSH`. <https://pubmed.ncbi.nlm.nih.gov/389116/>
64. Sinha P, et al. **Physiologic analysis and clinical performance of the ventilatory ratio in ARDS.** *Am J Respir Crit Care Med* 2019;199:333-41. <https://pubmed.ncbi.nlm.nih.gov/30211618/>

## 12. The right ventricle

65. 🔑 Boissier F, et al. **Prevalence and prognosis of cor pulmonale during protective ventilation for ARDS.** *Intensive Care Med* 2013;39:1725-33. — acute cor pulmonale in ~22%; the model's `AcuteCorPulm`. <https://pubmed.ncbi.nlm.nih.gov/23673401/>
66. Vieillard-Baron A, et al. **Acute cor pulmonale in ARDS submitted to protective ventilation: incidence, clinical implications and prognosis.** *Crit Care Med* 2001;29:1551-5. <https://pubmed.ncbi.nlm.nih.gov/11505125/>
67. Repessé X, Vieillard-Baron A. **Right heart function during acute respiratory distress syndrome.** *Ann Transl Med* 2017;5:295. <https://pubmed.ncbi.nlm.nih.gov/28828367/>

## 13. Fibroproliferation vs resolution

68. Marshall RP, et al. **Fibroproliferation occurs early in ARDS and impacts on outcome.** *Am J Respir Crit Care Med* 2000;162:1783-8. — N-terminal type III procollagen peptide; the model's `COLL`. <https://pubmed.ncbi.nlm.nih.gov/11069813/>
69. 🔑 Pittet JF, et al. **TGF-β is a critical mediator of acute lung injury.** *J Clin Invest* 2001;107:1537-44. — the αvβ6-TGF-β axis; the model's `TGFB`. <https://pubmed.ncbi.nlm.nih.gov/11413161/>
70. Burnham EL, et al. **The fibroproliferative response in ARDS: mechanisms and clinical significance.** *Eur Respir J* 2014;43:276-85. <https://pubmed.ncbi.nlm.nih.gov/23520315/>
71. Levy BD, Serhan CN. **Resolution of acute inflammation in the lung.** *Annu Rev Physiol* 2014;76:467-92. — lipoxins and resolvins; the model's M2/efferocytosis terms. <https://pubmed.ncbi.nlm.nih.gov/24313723/>
72. D'Alessio FR, et al. **CD4+CD25+Foxp3+ Tregs resolve experimental lung injury in mice.** *J Clin Invest* 2009;119:2898-913. <https://pubmed.ncbi.nlm.nih.gov/19770521/>

## 14. Subphenotypes and heterogeneity of treatment effect

73. 🔑 Calfee CS, et al. **Subphenotypes in acute respiratory distress syndrome: latent class analysis of data from two randomised controlled trials.** *Lancet Respir Med* 2014;2:611-20. — the hyperinflammatory/hypoinflammatory latent classes; the model's `HYPER` parameter. <https://pubmed.ncbi.nlm.nih.gov/24853585/>
74. 🔑 Famous KR, et al. **ARDS subphenotypes respond differently to randomized fluid management strategy.** *Am J Respir Crit Care Med* 2017;195:331-8. — opposing effects of the FACTT fluid strategy by subphenotype. <https://pubmed.ncbi.nlm.nih.gov/27513822/>
75. Calfee CS, et al. **Acute respiratory distress syndrome subphenotypes and differential response to simvastatin (HARP-2 secondary analysis).** *Lancet Respir Med* 2018;6:691-8. <https://pubmed.ncbi.nlm.nih.gov/30078618/>
76. Sinha P, et al. **Development and validation of parsimonious algorithms to classify ARDS phenotypes.** *Lancet Respir Med* 2020;8:247-57. — a three-variable classifier of IL-8, bicarbonate, and sTNFR1. <https://pubmed.ncbi.nlm.nih.gov/31948926/>
77. Constantin JM, et al. **Personalised mechanical ventilation tailored to lung morphology versus low PEEP in patients with ARDS (LIVE).** *Lancet Respir Med* 2019;7:870-80. — focal versus non-focal; harm when the strategy is mismatched. <https://pubmed.ncbi.nlm.nih.gov/31399381/>
78. Gattinoni L, et al. **COVID-19 pneumonia: different respiratory treatments for different phenotypes?** *Intensive Care Med* 2020;46:1099-102. — the L and H types. <https://pubmed.ncbi.nlm.nih.gov/32291463/>

## 15. Randomised trials of ventilation strategies

79. 🔑 **Acute Respiratory Distress Syndrome Network (ARMA). Ventilation with lower tidal volumes as compared with traditional tidal volumes for acute lung injury and ARDS.** *N Engl J Med* 2000;342:1301-8. — 39.8% → 31.0%. The primary target for calibrating the model's hazard. <https://pubmed.ncbi.nlm.nih.gov/10793162/>
80. **ARDSNet (ALVEOLI). Higher versus lower positive end-expiratory pressures in patients with ARDS.** *N Engl J Med* 2004;351:327-36. <https://pubmed.ncbi.nlm.nih.gov/15269312/>
81. Meade MO, et al. **Ventilation strategy using low tidal volumes, recruitment maneuvers, and high PEEP (LOVS).** *JAMA* 2008;299:637-45. <https://pubmed.ncbi.nlm.nih.gov/18270352/>
82. Mercat A, et al. **Positive end-expiratory pressure setting in adults with acute lung injury and ARDS (EXPRESS).** *JAMA* 2008;299:646-55. <https://pubmed.ncbi.nlm.nih.gov/18270353/>
83. 🔑 **Writing Group for the Alveolar Recruitment for ARDS Trial (ART) Investigators. Effect of lung recruitment and titrated PEEP vs low PEEP on mortality in patients with moderate to severe ARDS.** *JAMA* 2017;318:1335-45. — aggressive recruitment increased mortality. Reproduced as the harm of high PEEP in the model. <https://pubmed.ncbi.nlm.nih.gov/28973363/>
84. Briel M, et al. **Higher vs lower PEEP in patients with ARDS: systematic review and meta-analysis.** *JAMA* 2010;303:865-73. <https://pubmed.ncbi.nlm.nih.gov/20197533/>
85. 🔑 Guérin C, et al. **Prone positioning in severe acute respiratory distress syndrome (PROSEVA).** *N Engl J Med* 2013;368:2159-68. — 32.8% → 16.0%. Scenarios 4/P2 of the model. <https://pubmed.ncbi.nlm.nih.gov/23688302/>
86. Guérin C, et al. **A prospective international observational prevalence study on prone positioning of ARDS patients (APRONET).** *Intensive Care Med* 2018;44:22-37. <https://pubmed.ncbi.nlm.nih.gov/29218379/>
87. Ehrmann S, et al. **Awake prone positioning for COVID-19 acute hypoxaemic respiratory failure: a randomised, controlled, multinational, open-label meta-trial.** *Lancet Respir Med* 2021;9:1387-95. <https://pubmed.ncbi.nlm.nih.gov/34425070/>

## 16. Neuromuscular blockade and sedation

88. 🔑 Papazian L, et al. **Neuromuscular blockers in early acute respiratory distress syndrome (ACURASYS).** *N Engl J Med* 2010;363:1107-16. — 48 hours of cisatracurium. <https://pubmed.ncbi.nlm.nih.gov/20843245/>
89. 🔑 **National Heart, Lung, and Blood Institute PETAL Clinical Trials Network (ROSE). Early neuromuscular blockade in the acute respiratory distress syndrome.** *N Engl J Med* 2019;380:1997-2008. — neutral against a lightly sedated control arm. The `SEDATION` dependence of the model. <https://pubmed.ncbi.nlm.nih.gov/31112383/>
90. Slutsky AS. **Neuromuscular blocking agents in ARDS.** *N Engl J Med* 2010;363:1176-80 (editorial). <https://pubmed.ncbi.nlm.nih.gov/20843254/>

## 17. Corticosteroids

91. 🔑 Villar J, et al. **Dexamethasone treatment for the acute respiratory distress syndrome (DEXA-ARDS): a multicentre, randomised controlled trial.** *Lancet Respir Med* 2020;8:267-76. — 20 mg × 5 days → 10 mg × 5 days; VFD +4.8, mortality 36% → 21%. Scenario 6 of the model. <https://pubmed.ncbi.nlm.nih.gov/32043986/>
92. 🔑 Steinberg KP, et al. **Efficacy and safety of corticosteroids for persistent acute respiratory distress syndrome (LaSRS).** *N Engl J Med* 2006;354:1671-84. — increased mortality when started after day 14. Scenario 13 of the model. <https://pubmed.ncbi.nlm.nih.gov/16625008/>
93. Meduri GU, et al. **Methylprednisolone infusion in early severe ARDS: results of a randomized controlled trial.** *Chest* 2007;131:954-63. <https://pubmed.ncbi.nlm.nih.gov/17426195/>
94. 🔑 **RECOVERY Collaborative Group. Dexamethasone in hospitalized patients with Covid-19.** *N Engl J Med* 2021;384:693-704. <https://pubmed.ncbi.nlm.nih.gov/32678530/>
95. Chaudhuri D, et al. **Corticosteroids in COVID-19 and non-COVID-19 ARDS: a systematic review and meta-analysis.** *Intensive Care Med* 2021;47:521-37. <https://pubmed.ncbi.nlm.nih.gov/33876268/>
96. Guo RF, Ward PA. **Mediators and regulation of neutrophil accumulation in inflammatory responses in lung** (background on glucocorticoid transrepression). *Immunol Res* 2002;24:275-86. <https://pubmed.ncbi.nlm.nih.gov/11817326/>

## 18. Immunomodulators in COVID-ARDS and beyond

97. 🔑 **RECOVERY Collaborative Group. Tocilizumab in patients admitted to hospital with COVID-19.** *Lancet* 2021;397:1637-45. <https://pubmed.ncbi.nlm.nih.gov/33933206/>
98. **REMAP-CAP Investigators. Interleukin-6 receptor antagonists in critically ill patients with Covid-19.** *N Engl J Med* 2021;384:1491-502. <https://pubmed.ncbi.nlm.nih.gov/33631065/>
99. Kalil AC, et al. **Baricitinib plus remdesivir for hospitalized adults with Covid-19 (ACTT-2).** *N Engl J Med* 2021;384:795-807. <https://pubmed.ncbi.nlm.nih.gov/33306283/>
100. Marconi VC, et al. **Efficacy and safety of baricitinib for the treatment of hospitalised adults with COVID-19 (COV-BARRIER).** *Lancet Respir Med* 2021;9:1407-18. <https://pubmed.ncbi.nlm.nih.gov/34480861/>
101. Matthay MA, et al. **Treatment with allogeneic mesenchymal stromal cells for moderate to severe ARDS (START phase 2a).** *Lancet Respir Med* 2019;7:154-62. <https://pubmed.ncbi.nlm.nih.gov/30455077/>

## 19. Negative pharmacological trials — why the model must not "cure" ARDS with drugs

102. 🔑 Taylor RW, et al. **Low-dose inhaled nitric oxide in patients with acute lung injury: a randomized controlled trial.** *JAMA* 2004;291:1603-9. — oxygenation improves, mortality does not. The validation target of scenario 8 of the model. <https://pubmed.ncbi.nlm.nih.gov/15069048/>
103. Gebistorf F, et al. **Inhaled nitric oxide for acute respiratory distress syndrome in children and adults.** *Cochrane Database Syst Rev* 2016;(6):CD002787. <https://pubmed.ncbi.nlm.nih.gov/27347773/>
104. 🔑 Gao Smith F, et al. **Effect of intravenous β-2 agonist treatment on clinical outcomes in acute respiratory distress syndrome (BALTI-2).** *Lancet* 2012;379:229-35. — stopped early (harm). The warning attached to the model's `BETA2` switch. <https://pubmed.ncbi.nlm.nih.gov/22166903/>
105. **National Heart, Lung, and Blood Institute ARDS Clinical Trials Network (SAILS). Rosuvastatin for sepsis-associated acute respiratory distress syndrome.** *N Engl J Med* 2014;370:2191-200. <https://pubmed.ncbi.nlm.nih.gov/24835849/>
106. McAuley DF, et al. **Simvastatin in the acute respiratory distress syndrome (HARP-2).** *N Engl J Med* 2014;371:1695-703. <https://pubmed.ncbi.nlm.nih.gov/25268516/>
107. Fowler AA 3rd, et al. **Effect of vitamin C infusion on organ failure and biomarkers of inflammation and vascular injury in patients with sepsis and severe acute respiratory failure (CITRIS-ALI).** *JAMA* 2019;322:1261-70. <https://pubmed.ncbi.nlm.nih.gov/31573637/>
108. **ARDS Network. Ketoconazole for early treatment of acute lung injury and ARDS.** *JAMA* 2000;283:1995-2002. <https://pubmed.ncbi.nlm.nih.gov/10789667/>

## 20. Extracorporeal membrane oxygenation (ECMO)

109. 🔑 Combes A, et al. **Extracorporeal membrane oxygenation for severe acute respiratory distress syndrome (EOLIA).** *N Engl J Med* 2018;378:1965-75. — 35% vs 46%. Scenario 9 of the model. <https://pubmed.ncbi.nlm.nih.gov/29791822/>
110. Peek GJ, et al. **Efficacy and economic assessment of conventional ventilatory support versus extracorporeal membrane oxygenation for severe adult respiratory failure (CESAR).** *Lancet* 2009;374:1351-63. <https://pubmed.ncbi.nlm.nih.gov/19762075/>
111. Goligher EC, et al. **Extracorporeal membrane oxygenation for severe ARDS and posterior probability of mortality benefit in a post hoc Bayesian analysis of a randomized clinical trial.** *JAMA* 2018;320:2251-9. <https://pubmed.ncbi.nlm.nih.gov/30347031/>
112. Schmidt M, et al. **Mechanical ventilation management during ECMO for ARDS: an international multicenter prospective cohort.** *Am J Respir Crit Care Med* 2019;200:1002-12. <https://pubmed.ncbi.nlm.nih.gov/31144997/>

## 21. Long-term outcomes

113. Herridge MS, et al. **Functional disability 5 years after acute respiratory distress syndrome.** *N Engl J Med* 2011;364:1293-304. <https://pubmed.ncbi.nlm.nih.gov/21470008/>
114. Herridge MS, et al. **One-year outcomes in survivors of the acute respiratory distress syndrome.** *N Engl J Med* 2003;348:683-93. <https://pubmed.ncbi.nlm.nih.gov/12594312/>
115. Pandharipande PP, et al. **Long-term cognitive impairment after critical illness.** *N Engl J Med* 2013;369:1306-16. <https://pubmed.ncbi.nlm.nih.gov/24088092/>
116. 🔑 Puthucheary ZA, et al. **Acute skeletal muscle wasting in critical illness.** *JAMA* 2013;310:1591-600. — the model's `WEAK` compartment. <https://pubmed.ncbi.nlm.nih.gov/24108501/>
117. Schweickert WD, et al. **Early physical and occupational therapy in mechanically ventilated, critically ill patients: a randomised controlled trial.** *Lancet* 2009;373:1874-82. — the model's `MOBILISE`. <https://pubmed.ncbi.nlm.nih.gov/19446324/>

## 22. Methodology — trial endpoints and QSP

118. Yehya N, Harhay MO, Curley MAQ, et al. **Reappraisal of ventilator-free days in critical care research.** *Am J Respir Crit Care Med* 2019;200:828-36. — the definition of VFD and its survival weighting; the model's `VFD_expected`. <https://pubmed.ncbi.nlm.nih.gov/31034248/>
119. Baron KT, Gastonguay MR. **Simulation from ODE-based population PK/PD and systems pharmacology models in R with mrgsolve.** *J Pharmacokinet Pharmacodyn* 2015;42:S84-5 (mrgsolve methodology). <https://pubmed.ncbi.nlm.nih.gov/26373673/>
120. Ramanan M, et al. **A systems pharmacology perspective on critical-care trial design** — a related review: Vincent JL, et al. **Why have most trials of new treatments for sepsis and ARDS failed?** *Lancet Respir Med* 2022;10:e64-5. <https://pubmed.ncbi.nlm.nih.gov/35803300/>

---

### Summary of the key figures used directly in calibrating the model

| Item | Value | Source |
|---|---|---|
| 28-day mortality, Vt 12 → 6 mL/kg | 39.8% → 31.0% | ARMA (#79) |
| VFD, Vt 12 → 6 mL/kg | 10 → 12 days | ARMA (#79) |
| 28-day mortality, prone positioning (severe) | 32.8% → 16.0% | PROSEVA (#85) |
| 28/60-day mortality, dexamethasone | 36% → 21%, VFD +4.8 | DEXA-ARDS (#91) |
| Conservative fluids | VFD +2.5, no change in mortality | FACTT (#46) |
| 60-day mortality, VV-ECMO (severe) | 46% → 35% | EOLIA (#109) |
| iNO | oxygenation improved, mortality unchanged | Taylor 2004 (#102) |
| Prevalence of acute cor pulmonale | 20-25% | Boissier 2013 (#65) |
| Criterion for impaired AFC | <3%/h | Ware 2001 (#8) |
| Vd/Vt risk of death | independently raised above 0.6 | Nuckton 2002 (#62) |
| Mechanical power threshold | ~12 J/min | Gattinoni 2016 (#50) |
| Driving pressure threshold | ΔP > 15 cmH2O | Amato 2015 (#49) |
