# References — Delayed cerebral ischaemia after aneurysmal subarachnoid haemorrhage (aSAH-DCI)
# References — Aneurysmal Subarachnoid Haemorrhage & Delayed Cerebral Ischaemia

This list gives the evidence for the structure that `sah_qsp_model.dot` ·
`sah_mrgsolve_model.R` · `sah_shiny_app.R` claim. The annotation on each entry records
"which term of the model this paper justifies".

**Every PMID was looked up through the PubMed E-utilities (esearch/esummary) to confirm its
title, journal, and year, and each title is reproduced exactly as PubMed returned it.** Candidate
references whose lookup failed were left out of the list (no PMID has been guessed).

## The three structural claims of this model and their evidence

| Claim | Content | Primary evidence |
|------|------|----------|
| CLAIM 1 | Four consumers share a single arteriolar reserve — large-vessel spasm · microvascular tone · microthrombosis · perfusion pressure/SD | Rowland 2012 · Macdonald 2014 · Etminan 2011 |
| CLAIM 2 | The microcirculation is not a footnote to resistance; it erodes the oxygen extraction ceiling (CTH → OEFmax↓) | Jespersen & Østergaard 2012 · Østergaard 2013 |
| CLAIM 3 | The time window is not a calendar but the kinetics of haemoglobin handling (a transition chain plus an inducible sink) | Suzuki 1999 (HO-1) · Hugelshofer 2019 (haptoglobin) |

The single most important contrast: reading **CONSCIOUS-1** (65% reduction in vasospasm) and
**CONSCIOUS-2** (no change in outcome) side by side, at the same drug and the same dose.
This model decomposes that gap as 'redundancy first, harm channel second' and
computes it.

---

**89 references in total** (all PMID-verified).

## 1. Overview · Epidemiology · Definitions & Scales

1. **Subarachnoid haemorrhage.** *Lancet* 2007. [PMID 17258671](https://pubmed.ncbi.nlm.nih.gov/17258671/)  
   A classic review. The starting point for this model's clinical skeleton (grades, time axis, complications).

2. **Spontaneous subarachnoid haemorrhage.** *Lancet* 2017. [PMID 27637674](https://pubmed.ncbi.nlm.nih.gov/27637674/)  
   A modern review. The representative paper for the view that describes DCI separately from vasospasm.

3. **2023 Guideline for the Management of Patients With Aneurysmal Subarachnoid Hemorrhage: A Guideline From the American Heart Association/American Stroke Association.** *Stroke* 2023. [PMID 37212182](https://pubmed.ncbi.nlm.nih.gov/37212182/)  
   The 2023 AHA/ASA guideline. Evidence for the nimodipine recommendation, the uncertainty around induced hypertension, and target blood pressure and fluid management.

4. **Comment on the 2023 Guidelines for the Management of Patients With Aneurysmal Subarachnoid Hemorrhage.** *Stroke* 2023. [PMID 37581267](https://pubmed.ncbi.nlm.nih.gov/37581267/)  
   A commentary on the guideline above. Shows the dispute over how to interpret the strength of the recommendations.

5. **Critical care management of patients following aneurysmal subarachnoid hemorrhage: recommendations from the Neurocritical Care Society's Multidisciplinary Consensus Conference.** *Neurocrit Care* 2011. [PMID 21773873](https://pubmed.ncbi.nlm.nih.gov/21773873/)  
   The Neurocritical Care Society consensus. Evidence for the use of monitoring (TCD·PbtO2·microdialysis).

6. **Definition of delayed cerebral ischemia after aneurysmal subarachnoid hemorrhage as an outcome event in clinical trials and observational studies: proposal of a multidisciplinary research group.** *Stroke* 2010. [PMID 20798370](https://pubmed.ncbi.nlm.nih.gov/20798370/)  
   The standard definition of DCI (neurological deterioration + new infarction). The model's DCI decision rule is an operationalisation of this definition.

7. **Report of World Federation of Neurological Surgeons Committee on a Universal Subarachnoid Hemorrhage Grading Scale.** *J Neurosurg* 1988. [PMID 3131498](https://pubmed.ncbi.nlm.nih.gov/3131498/)  
   The original WFNS grading paper. The model's initial value for EBIB (early brain injury burden) is derived from WFNS.

8. **Prediction of symptomatic vasospasm after subarachnoid hemorrhage: the modified fisher scale.** *Neurosurgery* 2006. [PMID 16823296](https://pubmed.ncbi.nlm.nih.gov/16823296/)  
   The modified Fisher scale. The model's CLOT0·IVH0 initial conditions are set from this grade.

9. **Reliability of the modified Rankin Scale applied by telephone.** *Neurol Int* 2013. [PMID 23717781](https://pubmed.ncbi.nlm.nih.gov/23717781/)  
   Reliability of the mRS. Background on measurement error for the endpoint definition of the outcome logit.


## 2. Vasospasm–DCI dissociation: the phenomenon this model sets out to explain

10. **Delayed neurological deterioration after subarachnoid haemorrhage.** *Nat Rev Neurol* 2014. [PMID 24323051](https://pubmed.ncbi.nlm.nih.gov/24323051/)  
   The key review arguing that delayed neurological deterioration cannot be explained by vasospasm alone. The problem this model starts from.

11. **Delayed cerebral ischaemia after subarachnoid haemorrhage: looking beyond vasospasm.** *Br J Anaesth* 2012. [PMID 22879655](https://pubmed.ncbi.nlm.nih.gov/22879655/)  
   Exactly as its title says, 'DCI: beyond vasospasm'. The descriptive precursor of the four-consumer structure.

12. **Relationship between vasospasm, cerebral perfusion, and delayed cerebral ischemia after aneurysmal subarachnoid hemorrhage.** *Neuroradiology* 2009. [PMID 19623472](https://pubmed.ncbi.nlm.nih.gov/19623472/)  
   A study showing quantitatively that the vasospasm–perfusion–DCI triangle is loose.

13. **Relationship between angiographic vasospasm and regional hypoperfusion in aneurysmal subarachnoid hemorrhage.** *Stroke* 2012. [PMID 22492520](https://pubmed.ncbi.nlm.nih.gov/22492520/)  
   A study showing that the relation between vasospasm and regional hypoperfusion is only partial. The limit of the large-vessel pathway's explanatory power.

14. **The relationship between delayed infarcts and angiographic vasospasm after aneurysmal subarachnoid hemorrhage.** *Neurosurgery* 2013. [PMID 23313984](https://pubmed.ncbi.nlm.nih.gov/23313984/)  
   Directly measures the mismatch between delayed infarction and angiographic spasm. Corresponds to the model's r² calculation.

15. **Angiographic vasospasm is strongly correlated with cerebral infarction after subarachnoid hemorrhage.** *Stroke* 2011. [PMID 21350201](https://pubmed.ncbi.nlm.nih.gov/21350201/)  
   Evidence in the opposite direction: reports a strong correlation between spasm and infarction. The model has to explain both.

16. **Predictors of cerebral infarction in aneurysmal subarachnoid hemorrhage.** *Stroke* 2004. [PMID 15218156](https://pubmed.ncbi.nlm.nih.gov/15218156/)  
   Predictors of cerebral infarction. The point that clinical grade and haemorrhage volume predict more strongly than spasm does.

17. **Effect of pharmaceutical treatment on vasospasm, delayed cerebral ischemia, and clinical outcome in patients with aneurysmal subarachnoid hemorrhage: a systematic review and meta-analysis.** *J Cereb Blood Flow Metab* 2011. [PMID 21285966](https://pubmed.ncbi.nlm.nih.gov/21285966/)  
   The decisive meta-analysis: across the whole body of trials, pharmacological treatment reduces vasospasm but fails to change outcome. The phenomenon this model sets out to explain, in itself.

18. **Effect of endothelin-receptor antagonists on delayed cerebral ischemia after aneurysmal subarachnoid hemorrhage remains unclear.** *Stroke* 2009. [PMID 19875735](https://pubmed.ncbi.nlm.nih.gov/19875735/)  
   Meta-analysis of endothelin receptor antagonists. A quantitative summary of reduced spasm vs unchanged outcome.


## 3. Clazosentan and the endothelin axis: an instructive failure

19. **Clazosentan to overcome neurological ischemia and infarction occurring after subarachnoid hemorrhage (CONSCIOUS-1): randomized, double-blind, placebo-controlled phase 2 dose-finding trial.** *Stroke* 2008. [PMID 18688013](https://pubmed.ncbi.nlm.nih.gov/18688013/)  
   CONSCIOUS-1. Clazosentan reduced moderate–severe vasospasm by about 65% (model calibration target RR 0.35).

20. **Clazosentan, an endothelin receptor antagonist, in patients with aneurysmal subarachnoid haemorrhage undergoing surgical clipping: a randomised, double-blind, placebo-controlled phase 3 trial (CONSCIOUS-2).** *Lancet Neurol* 2011. [PMID 21640651](https://pubmed.ncbi.nlm.nih.gov/21640651/)  
   CONSCIOUS-2. The same drug failed to improve outcome at 3 months (model calibration target RR ~1.05). The central case for this model.

21. **Randomized trial of clazosentan in patients with aneurysmal subarachnoid hemorrhage undergoing endovascular coiling.** *Stroke* 2012. [PMID 22403047](https://pubmed.ncbi.nlm.nih.gov/22403047/)  
   CONSCIOUS-3. Confirmation in coiled patients. At high dose the harm channel — pulmonary oedema, hypotension — is pronounced.

22. **Effects of clazosentan on cerebral vasospasm-related morbidity and all-cause mortality after aneurysmal subarachnoid hemorrhage: two randomized phase 3 trials in Japanese patients.** *J Neurosurg* 2022. [PMID 35364589](https://pubmed.ncbi.nlm.nih.gov/35364589/)  
   The Japanese phase 3 trial. It reduced vasospasm-related morbidity and mortality, showing that the conclusion changes with the endpoint definition.

23. **The REACT study: design of a randomized phase 3 trial to assess the efficacy and safety of clazosentan for preventing deterioration due to delayed cerebral ischemia after aneurysmal subarachnoid hemorrhage.** *BMC Neurol* 2022. [PMID 36539711](https://pubmed.ncbi.nlm.nih.gov/36539711/)  
   The REACT trial design. Where the re-evaluation of clazosentan now stands.

24. **Endothelin and subarachnoid hemorrhage: an overview.** *Neurosurgery* 1998. [PMID 9766314](https://pubmed.ncbi.nlm.nih.gov/9766314/)  
   A review of the endothelin axis. The ET-1/ETA/ETB structure and the model's ETeff term.

25. **Plasma endothelin-1 as screening marker for cerebral vasospasm after subarachnoid hemorrhage.** *Neurocrit Care* 2014. [PMID 23921571](https://pubmed.ncbi.nlm.nih.gov/23921571/)  
   The screening performance of plasma ET-1. The loose link between effector concentration and clinical events.


## 4. Nimodipine and calcium antagonists

26. **Effect of oral nimodipine on cerebral infarction and outcome after subarachnoid haemorrhage: British aneurysm nimodipine trial.** *BMJ* 1989. [PMID 2496789](https://pubmed.ncbi.nlm.nih.gov/2496789/)  
   BRANT. Oral nimodipine reduced cerebral infarction and poor outcome. The only drug that changed outcome.

27. **Cerebral arterial spasm--a controlled trial of nimodipine in patients with subarachnoid hemorrhage.** *N Engl J Med* 1983. [PMID 6338383](https://pubmed.ncbi.nlm.nih.gov/6338383/)  
   The first nimodipine RCT. The prototype case in which angiographic spasm barely changed yet clinical outcome improved.

28. **Calcium antagonists for aneurysmal subarachnoid haemorrhage.** *Cochrane Database Syst Rev* 2007. [PMID 17636626](https://pubmed.ncbi.nlm.nih.gov/17636626/)  
   The Cochrane review of calcium antagonists. RR for poor outcome about 0.67 (model calibration target).

29. **Nimodipine-Induced Blood Pressure Changes Can Predict Delayed Cerebral Ischemia.** *Front Neurol* 2019. [PMID 31736865](https://pubmed.ncbi.nlm.nih.gov/31736865/)  
   The blood pressure change caused by nimodipine itself predicts DCI. Direct evidence for the model's CPP harm channel (DMAP_NIM).

30. **NEWTON: Nimodipine Microparticles to Enhance Recovery While Reducing Toxicity After Subarachnoid Hemorrhage.** *Neurocrit Care* 2015. [PMID 25678453](https://pubmed.ncbi.nlm.nih.gov/25678453/)  
   Phase 1/2 of intraventricular sustained-release nimodipine (EG-1962). The strategy of raising local exposure without systemic hypotension.

31. **NEWTON-2 Cisternal (Nimodipine Microparticles to Enhance Recovery While Reducing Toxicity After Subarachnoid Hemorrhage): A Phase 2, Multicenter, Randomized, Open-Label Safety Study of Intracisternal EG-1962 in Aneurysmal Subarachnoid Hemorrhage.** *Neurosurgery* 2020. [PMID 32985652](https://pubmed.ncbi.nlm.nih.gov/32985652/)  
   NEWTON-2. Another case in which local delivery reduced vasospasm without improving outcome.

32. **Nicardipine prolonged-release implants for preventing cerebral vasospasm after subarachnoid hemorrhage: effect and outcome in the first 100 patients.** *Neurol Med Chir (Tokyo)* 2007. [PMID 17895611](https://pubmed.ncbi.nlm.nih.gov/17895611/)  
   Nicardipine sustained-release implants. A strong effect of local large-vessel relaxation.

33. **Intrathecal nicardipine for cerebral vasospasm after non-traumatic subarachnoid hemorrhage: a meta-analysis.** *Neurosurg Rev* 2025. [PMID 40295406](https://pubmed.ncbi.nlm.nih.gov/40295406/)  
   Meta-analysis of intrathecal nicardipine. Evidence for the model's S9 (intrathecal nicardipine) scenario.


## 5. Other interventional trials

34. **Simvastatin in aneurysmal subarachnoid haemorrhage (STASH): a multicentre randomised phase 3 trial.** *Lancet Neurol* 2014. [PMID 24837690](https://pubmed.ncbi.nlm.nih.gov/24837690/)  
   STASH. Simvastatin failed to improve outcome (model calibration target RR ~1.0).

35. **Magnesium for aneurysmal subarachnoid haemorrhage (MASH-2): a randomised placebo-controlled trial.** *Lancet* 2012. [PMID 22633825](https://pubmed.ncbi.nlm.nih.gov/22633825/)  
   MASH-2. Magnesium ineffective.

36. **Intravenous magnesium sulphate for aneurysmal subarachnoid hemorrhage (IMASH): a randomized, double-blinded, placebo-controlled, multicenter phase III trial.** *Stroke* 2010. [PMID 20378868](https://pubmed.ncbi.nlm.nih.gov/20378868/)  
   IMASH. Magnesium ineffective (confirmation).

37. **Effects of cilostazol on cerebral vasospasm after aneurysmal subarachnoid hemorrhage: a multicenter prospective, randomized, open-label blinded end point trial.** *J Neurosurg* 2013. [PMID 23039152](https://pubmed.ncbi.nlm.nih.gov/23039152/)  
   A cilostazol RCT. A multi-pathway action: antiplatelet plus microvascular dilatation.

38. **Efficacy of Cilostazol in Prevention of Delayed Cerebral Ischemia after Aneurysmal Subarachnoid Hemorrhage: A Meta-Analysis.** *J Stroke Cerebrovasc Dis* 2018. [PMID 30093204](https://pubmed.ncbi.nlm.nih.gov/30093204/)  
   Meta-analysis of cilostazol. RR for DCI about 0.47 (model calibration target).

39. **Effectiveness of Lumbar Cerebrospinal Fluid Drain Among Patients With Aneurysmal Subarachnoid Hemorrhage: A Randomized Clinical Trial.** *JAMA Neurol* 2023. [PMID 37330974](https://pubmed.ncbi.nlm.nih.gov/37330974/)  
   EARLYDRAIN. Early lumbar drainage reduced poor outcome. The only intervention in the model that 'truncates the input function'.

40. **Induced Hypertension for Delayed Cerebral Ischemia After Aneurysmal Subarachnoid Hemorrhage: A Randomized Clinical Trial.** *Stroke* 2018. [PMID 29158449](https://pubmed.ncbi.nlm.nih.gov/29158449/)  
   HIMALAIA. The randomised trial of induced hypertension was ineffective. The model explains this as an effect that splits by autoregulatory state.

41. **Ultra-early tranexamic acid after subarachnoid haemorrhage (ULTRA): a randomised controlled trial.** *Lancet* 2021. [PMID 33357465](https://pubmed.ncbi.nlm.nih.gov/33357465/)  
   ULTRA. Ultra-early tranexamic acid reduced rebleeding but did not improve outcome.

42. **Milrinone and homeostasis to treat cerebral vasospasm associated with subarachnoid hemorrhage: the Montreal Neurological Hospital protocol.** *Neurocrit Care* 2012. [PMID 22528278](https://pubmed.ncbi.nlm.nih.gov/22528278/)  
   The Montreal milrinone protocol. Rescue therapy combining microvascular dilatation with increased cardiac output.

43. **Effect of prophylactic transluminal balloon angioplasty on cerebral vasospasm and outcome in patients with Fisher grade III subarachnoid hemorrhage: results of a phase II multicenter, randomized, clinical trial.** *Stroke* 2008. [PMID 18420953](https://pubmed.ncbi.nlm.nih.gov/18420953/)  
   Prophylactic balloon angioplasty. A case in which physically widening the large vessels did not improve outcome.

44. **Endovascular Rescue Treatment for Delayed Cerebral Ischemia After Subarachnoid Hemorrhage Is Safe and Effective.** *Front Neurol* 2019. [PMID 30858818](https://pubmed.ncbi.nlm.nih.gov/30858818/)  
   Endovascular rescue therapy. The last resort for intervening on the large-vessel pathway in real practice.


## 6. Haemoglobin · haptoglobin · HO-1: the clock

45. **A review of hemoglobin and the pathogenesis of cerebral vasospasm.** *Stroke* 1991. [PMID 1866764](https://pubmed.ncbi.nlm.nih.gov/1866764/)  
   A review of haemoglobin and the pathogenesis of vasospasm. Evidence for the model's OXYHB-centred structure.

46. **Haptoglobin administration into the subarachnoid space prevents hemoglobin-induced cerebral vasospasm.** *J Clin Invest* 2019. [PMID 31454333](https://pubmed.ncbi.nlm.nih.gov/31454333/)  
   Administering haptoglobin into the subarachnoid space prevents vasospasm. Demonstrates that scavenging capacity is the rate-limiting step.

47. **Haptoglobin genotype and outcome after aneurysmal subarachnoid haemorrhage.** *J Neurol Neurosurg Psychiatry* 2020. [PMID 31937585](https://pubmed.ncbi.nlm.nih.gov/31937585/)  
   Haptoglobin genotype and outcome. The model's HP22 covariate.

48. **Haptoglobin phenotype predicts cerebral vasospasm and clinical deterioration after aneurysmal subarachnoid hemorrhage.** *J Stroke Cerebrovasc Dis* 2013. [PMID 23498376](https://pubmed.ncbi.nlm.nih.gov/23498376/)  
   Haptoglobin phenotype predicts vasospasm and clinical deterioration.

49. **Heme oxygenase-1 gene induction as an intrinsic regulation against delayed cerebral vasospasm in rats.** *J Clin Invest* 1999. [PMID 10393699](https://pubmed.ncbi.nlm.nih.gov/10393699/)  
   HO-1 induction is the endogenous defence against delayed vasospasm. The model's 'inducible sink' and its delay (t½ 1.9 d).

50. **Bilirubin oxidation products (BOXes) and their role in cerebral vasospasm after subarachnoid hemorrhage.** *J Cereb Blood Flow Metab* 2006. [PMID 16467784](https://pubmed.ncbi.nlm.nih.gov/16467784/)  
   The vasoconstrictor action of bilirubin oxidation products (BOXes). The model's BOX state.

51. **Iron and early brain injury after subarachnoid hemorrhage.** *J Cereb Blood Flow Metab* 2010. [PMID 20736954](https://pubmed.ncbi.nlm.nih.gov/20736954/)  
   Iron and early brain injury. Evidence for the Fenton chemistry and ROS terms.

52. **Ferroptosis mechanisms in early brain injury after subarachnoid hemorrhage.** *Brain Res* 2026. [PMID 41802694](https://pubmed.ncbi.nlm.nih.gov/41802694/)  
   Haem/iron-mediated ferroptosis. The most recent synthesis of the oxidative injury pathway.


## 7. NO · endothelin · Rho-kinase: effector signalling

53. **Acute decrease in cerebral nitric oxide levels after subarachnoid hemorrhage.** *J Cereb Blood Flow Metab* 2000. [PMID 10724124](https://pubmed.ncbi.nlm.nih.gov/10724124/)  
   The acute fall in cerebral NO concentration immediately after SAH. The model's NOB term.

54. **Delayed cerebral vasospasm and nitric oxide: review, new hypothesis, and proposed treatment.** *Pharmacol Ther* 2005. [PMID 15626454](https://pubmed.ncbi.nlm.nih.gov/15626454/)  
   The NO deficiency hypothesis and a therapeutic proposal. NO scavenging by oxyHb is central.

55. **Nitric oxide in early brain injury after subarachnoid hemorrhage.** *Acta Neurochir Suppl* 2011. [PMID 21116923](https://pubmed.ncbi.nlm.nih.gov/21116923/)  
   NO in early brain injury. The link between EBI and vascular dysfunction.

56. **Involvement of Rho-kinase-mediated phosphorylation of myosin light chain in enhancement of cerebral vasospasm.** *Circ Res* 2000. [PMID 10926869](https://pubmed.ncbi.nlm.nih.gov/10926869/)  
   Rho-kinase-mediated MLC phosphorylation amplifies vasospasm. The model's RHOK (Ca sensitisation) term — why calcium blockade alone is not enough.

57. **Effect of AT877 on cerebral vasospasm after aneurysmal subarachnoid hemorrhage. Results of a prospective placebo-controlled double-blind trial.** *J Neurosurg* 1992. [PMID 1545249](https://pubmed.ncbi.nlm.nih.gov/1545249/)  
   The randomised trial of fasudil (AT877). Clinical evidence for Rho-kinase inhibition.


## 8. Microcirculation · transit time heterogeneity (CTH) · microthrombosis

58. **The role of the microcirculation in delayed cerebral ischemia and chronic degenerative changes after subarachnoid hemorrhage.** *J Cereb Blood Flow Metab* 2013. [PMID 24064495](https://pubmed.ncbi.nlm.nih.gov/24064495/)  
   The role the microcirculation plays in DCI. The primary source for CLAIM 2 of this model.

59. **The roles of cerebral blood flow, capillary transit time heterogeneity, and oxygen tension in brain oxygenation and metabolism.** *J Cereb Blood Flow Metab* 2012. [PMID 22044867](https://pubmed.ncbi.nlm.nih.gov/22044867/)  
   The relation between CBF · transit time heterogeneity (CTH) · oxygen tension. The model's form OEFmax = 0.85/(1+K·CTH) comes from here.

60. **Are We Barking Up the Wrong Vessels? Cerebral Microcirculation After Subarachnoid Hemorrhage.** *Stroke* 2015. [PMID 26152299](https://pubmed.ncbi.nlm.nih.gov/26152299/)  
   'Are we looking at the wrong vessels?' A manifesto review from the microcirculatory point of view.

61. **Capillary flow disturbances after experimental subarachnoid hemorrhage: A contributor to delayed cerebral ischemia?.** *Microcirculation* 2019. [PMID 30431201](https://pubmed.ncbi.nlm.nih.gov/30431201/)  
   Capillary flow disturbance after experimental SAH contributes to DCI. The experimental evidence for CTH.

62. **Microthrombosis after aneurysmal subarachnoid hemorrhage: an additional explanation for delayed cerebral ischemia.** *J Cereb Blood Flow Metab* 2008. [PMID 18628782](https://pubmed.ncbi.nlm.nih.gov/18628782/)  
   The representative paper arguing that microthrombosis is an additional explanation for DCI. The model's MTHR state.

63. **Mechanisms of microthrombosis and microcirculatory constriction after experimental subarachnoid hemorrhage.** *Acta Neurochir Suppl* 2013. [PMID 22890667](https://pubmed.ncbi.nlm.nih.gov/22890667/)  
   The mechanism of microthrombus formation and microcirculatory constriction.

64. **Intraoperative detection of early microvasospasm in patients with subarachnoid hemorrhage by using orthogonal polarization spectral imaging.** *Neurosurgery* 2003. [PMID 12762876](https://pubmed.ncbi.nlm.nih.gov/12762876/)  
   Direct intraoperative observation of early microvascular spasm. A pathway independent of large-vessel spasm.

65. **Reduction of neutrophil activity decreases early microvascular injury after subarachnoid haemorrhage.** *J Neuroinflammation* 2011. [PMID 21854561](https://pubmed.ncbi.nlm.nih.gov/21854561/)  
   Reducing neutrophil activity lessens early microvascular injury. The inflammation–microcirculation link.

66. **Pericyte constriction after stroke: the jury is still out.** *Nat Med* 2010. [PMID 20823870](https://pubmed.ncbi.nlm.nih.gov/20823870/)  
   The debate over pericyte constriction after stroke.

67. **Capillary pericytes regulate cerebral blood flow in health and disease.** *Nature* 2014. [PMID 24670647](https://pubmed.ncbi.nlm.nih.gov/24670647/)  
   The Nature study showing that capillary pericytes regulate CBF. The cellular basis of microvascular resistance.


## 9. Spreading depolarisation

68. **Delayed ischaemic neurological deficits after subarachnoid haemorrhage are associated with clusters of spreading depolarizations.** *Brain* 2006. [PMID 17067993](https://pubmed.ncbi.nlm.nih.gov/17067993/)  
   Delayed ischaemic deficits after SAH are associated with clusters of spreading depolarisation. The primary source for the model's SD pathway.

69. **The role of spreading depression, spreading depolarization and spreading ischemia in neurological disease.** *Nat Med* 2011. [PMID 21475241](https://pubmed.ncbi.nlm.nih.gov/21475241/)  
   A review of spreading depolarisation / spreading ischaemia. The conceptual basis for inverse neurovascular coupling.

70. **Recording, analysis, and interpretation of spreading depolarizations in neurointensive care: Review and recommendations of the COSBID research group.** *J Cereb Blood Flow Metab* 2017. [PMID 27317657](https://pubmed.ncbi.nlm.nih.gov/27317657/)  
   The COSBID recommendations. The standard for recording and interpreting SD.

71. **Delayed cerebral ischemia and spreading depolarization in absence of angiographic vasospasm after subarachnoid hemorrhage.** *J Cereb Blood Flow Metab* 2012. [PMID 22146193](https://pubmed.ncbi.nlm.nih.gov/22146193/)  
   DCI and SD in the absence of angiographic spasm. The direct evidence for the model treating SD as an independent consumer.

72. **Spreading Depolarizations: A Therapeutic Target Against Delayed Cerebral Ischemia After Subarachnoid Hemorrhage.** *J Clin Neurophysiol* 2016. [PMID 27258442](https://pubmed.ncbi.nlm.nih.gov/27258442/)  
   The view of SD as a therapeutic target. Evidence for the ketamine scenario.

73. **Preliminary evidence that ketamine inhibits spreading depolarizations in acute human brain injury.** *Stroke* 2009. [PMID 19520992](https://pubmed.ncbi.nlm.nih.gov/19520992/)  
   Preliminary evidence that ketamine suppresses SD in the human brain. The model's S12 scenario.


## 10. Cerebral autoregulation

74. **Clinical relevance of cerebral autoregulation following subarachnoid haemorrhage.** *Nat Rev Neurol* 2013. [PMID 23419369](https://pubmed.ncbi.nlm.nih.gov/23419369/)  
   The clinical significance of autoregulation after SAH. The model's AREG state and pressure passivity.

75. **Impairment of cerebral autoregulation predicts delayed cerebral ischemia after subarachnoid hemorrhage: a prospective observational study.** *Stroke* 2012. [PMID 23150652](https://pubmed.ncbi.nlm.nih.gov/23150652/)  
   Impaired autoregulation predicts DCI. The evidence for treating AREG as a risk factor, and the counterpart of the stratification result in R5.


## 11. Monitoring and biomarkers

76. **Transcranial Doppler versus angiography in patients with vasospasm due to a ruptured cerebral aneurysm: A systematic review.** *Stroke* 2001. [PMID 11588316](https://pubmed.ncbi.nlm.nih.gov/11588316/)  
   A systematic review of TCD vs angiography. The comparator for the model's calculation of TCD diagnostic performance.

77. **How is vasospasm screening using transcranial Doppler associated with delayed cerebral ischemia and outcomes in aneurysmal subarachnoid hemorrhage?.** *Acta Neurochir (Wien)* 2019. [PMID 30637487](https://pubmed.ncbi.nlm.nih.gov/30637487/)  
   The association of TCD screening with DCI and outcome. The point that flow velocity is confounded by flow.

78. **Monitoring of brain tissue oxygenation following severe subarachnoid hemorrhage.** *Neurol Res* 2003. [PMID 12866190](https://pubmed.ncbi.nlm.nih.gov/12866190/)  
   Brain tissue oxygen monitoring in severe SAH. Why PbtO2 measures the shared node directly.

79. **Cerebral metabolism and intracranial hypertension in high grade aneurysmal subarachnoid haemorrhage patients.** *Acta Neurochir Suppl* 2005. [PMID 16463827](https://pubmed.ncbi.nlm.nih.gov/16463827/)  
   Microdialysis metabolic markers and intracranial pressure. The evidence for a rising LPR.


## 12. Early brain injury and systemic complications

80. **Metamorphosis of subarachnoid hemorrhage research: from delayed vasospasm to early brain injury.** *Mol Neurobiol* 2011. [PMID 21161614](https://pubmed.ncbi.nlm.nih.gov/21161614/)  
   The research paradigm has moved from delayed vasospasm to early brain injury. Why EBI is the largest determinant of outcome in the model.

81. **Early Brain Injury After Subarachnoid Hemorrhage: Incidence and Mechanisms.** *Stroke* 2023. [PMID 36866673](https://pubmed.ncbi.nlm.nih.gov/36866673/)  
   The incidence and mechanisms of early brain injury (a 2023 synthesis).

82. **The incidence and pathophysiology of hyponatraemia after subarachnoid haemorrhage.** *Clin Endocrinol (Oxf)* 2006. [PMID 16487432](https://pubmed.ncbi.nlm.nih.gov/16487432/)  
   The frequency and mechanism of hyponatraemia after SAH. Calibration of the model's NAS state (<135, about 30-50%).

83. **Higher hemoglobin is associated with improved outcome after subarachnoid hemorrhage.** *Crit Care Med* 2007. [PMID 17717494](https://pubmed.ncbi.nlm.nih.gov/17717494/)  
   Higher haemoglobin is associated with better outcome. The clinical evidence that CaO2 enters the same node.

84. **Red blood cell transfusion in patients with subarachnoid hemorrhage: a multidisciplinary North American survey.** *Crit Care* 2011. [PMID 21244675](https://pubmed.ncbi.nlm.nih.gov/21244675/)  
   Red cell transfusion practice in SAH. The clinical context of the S13 scenario.

85. **Impact of cardiac complications on outcome after aneurysmal subarachnoid hemorrhage: a meta-analysis.** *Neurology* 2009. [PMID 19221297](https://pubmed.ncbi.nlm.nih.gov/19221297/)  
   The effect of cardiac complications on outcome (a meta-analysis). The neurogenic myocardial injury channel.

86. **Fever after subarachnoid hemorrhage: risk factors and impact on outcome.** *Neurology* 2007. [PMID 17314332](https://pubmed.ncbi.nlm.nih.gov/17314332/)  
   The effect of fever on outcome. In the model fever is a consumer on the CMRO2 (demand) side.

87. **Risk factors for surgery-related cerebral infarction after aneurysmal subarachnoid hemorrhage.** *Neurosurg Rev* 2025. [PMID 40316796](https://pubmed.ncbi.nlm.nih.gov/40316796/)  
   Risk factors for procedure-related cerebral infarction. The non-DCI component of infarct burden.


## 13. QSP methodology

88. **Quantitative Systems Pharmacology and Physiologically-Based Pharmacokinetic Modeling With mrgsolve: A Hands-On Tutorial.** *CPT Pharmacometrics Syst Pharmacol* 2019. [PMID 31652028](https://pubmed.ncbi.nlm.nih.gov/31652028/)  
   A hands-on QSP/PBPK tutorial using mrgsolve. The implementation tool for this model.

89. **Quantitative Systems Pharmacology: A Case for Disease Models.** *Clin Pharmacol Ther* 2017. [PMID 27709613](https://pubmed.ncbi.nlm.nih.gov/27709613/)  
   A discussion of QSP as disease modelling. The methodological position of this library.



---

## Calibration targets and model output values

All the values below were computed by `sah_reference_check.py` (an independent numpy/RK4 implementation), and
the full output is committed verbatim in `sah_reference_check_output.txt`.

| Item | Literature target | Source |
|------|-----------|------|
| Moderate–severe angiographic vasospasm (control arm) | about 66% | CONSCIOUS-1 [PMID 18688013](https://pubmed.ncbi.nlm.nih.gov/18688013/) |
| DCI incidence (standard care) | about 30% | Vergouwen 2010 [PMID 20798370](https://pubmed.ncbi.nlm.nih.gov/20798370/) |
| Poor outcome mRS 4-6 (90 days) | 30–35% | Macdonald 2017 [PMID 27637674](https://pubmed.ncbi.nlm.nih.gov/27637674/) |
| Hyponatraemia < 135 mmol/L | 30–50% | Sherlock 2006 [PMID 16487432](https://pubmed.ncbi.nlm.nih.gov/16487432/) |
| Clazosentan RR, moderate–severe spasm | 0.35 | CONSCIOUS-1 [PMID 18688013](https://pubmed.ncbi.nlm.nih.gov/18688013/) |
| Clazosentan RR, poor outcome | about 1.05 | CONSCIOUS-2 [PMID 21640651](https://pubmed.ncbi.nlm.nih.gov/21640651/) |
| Nimodipine RR, poor outcome | 0.67 | Cochrane [PMID 17636626](https://pubmed.ncbi.nlm.nih.gov/17636626/) |
| Cilostazol RR, DCI | about 0.47 | meta-analysis [PMID 30093204](https://pubmed.ncbi.nlm.nih.gov/30093204/) |
| Simvastatin RR, poor outcome | about 1.00 | STASH [PMID 24837690](https://pubmed.ncbi.nlm.nih.gov/24837690/) |
| Lumbar drainage RR, poor outcome | 0.76 | EARLYDRAIN [PMID 37330974](https://pubmed.ncbi.nlm.nih.gov/37330974/) |

## Where the model is shallower than, or differs from, the literature (recorded mismatches)

These are items left explicitly in place rather than tuned away. The full list and discussion are in
the "Where the model diverges from the literature or fails" section of [README.md](README.md) in the model's own directory.

1. **The model's outcome improvement with nimodipine is shallow.** Model RR (poor outcome) 0.80 vs Cochrane 0.67.
2. **The proportion of DCI without angiographic spasm is low.** Model 10.9% vs 20-30% reported in the literature.
3. **The two extremes of the modified Fisher scale are spread too far apart.** mF1 2.3% (literature about 6-12%),
   mF4 50.0% (literature about 35-40%). The ordering is correctly monotonic but the gradient is steep.
4. **Hyponatraemia < 130 mmol/L is low.** Model 4.6% vs 10-20% in the literature.
5. **A hypothesis that was not reproduced.** The prediction that "induced hypertension works only in pressure-passive patients"
   came out the other way round in the calculation (DCI RR 0.73 in the impaired-autoregulation group vs 0.62 in the intact group).
   The hypothesis that the HIMALAIA null result can be explained by autoregulatory stratification is not supported in this model.

## Methodology tools

- mrgsolve (R): <https://mrgsolve.org>
- gPKPDviz — a Shiny tool for mrgsolve-based PK/PD simulation
  - Paper: <https://pmc.ncbi.nlm.nih.gov/articles/PMC10941578/>
  - Code: <https://github.com/Genentech/gPKPDviz/>
- PubMed E-utilities (used to verify the PMIDs in this list): <https://www.ncbi.nlm.nih.gov/books/NBK25501/>
