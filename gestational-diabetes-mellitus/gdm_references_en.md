# Gestational Diabetes Mellitus (GDM) QSP Model — References
# Gestational Diabetes Mellitus — Annotated Reference List

> **PMID verification method.** Every PMID in this list was confirmed with NCBI E-utilities
> (`ecitmatch` citation matching + `esummary`) by checking journal, year, volume, pages, and first author
> against the record. The titles and bibliographic details are taken verbatim from the PubMed
> record, and citations that could not be confirmed were not put in this list.
>
> **All PMIDs below were machine-verified** against PubMed via E-utilities
> (`ecitmatch` citation matching, then `esummary` to confirm journal, year,
> volume, pages and first author). Titles are taken verbatim from the PubMed
> record. Citations that failed verification were dropped rather than guessed.

**85 references** in total · 12 sections. The square brackets at the end of each entry indicate which part of the model
that reference calibrates (for example `[→ CID weighting]`).

---

## 1. Landmark trials and large observational cohorts

Almost all of the clinical endpoint calibration of the model comes out of this section.

1. HAPO Study Cooperative Research Group. **Hyperglycemia and adverse pregnancy outcomes.** N Engl J Med 2008;358:1991-2002.
   <https://pubmed.ncbi.nlm.nih.gov/18463375/>
   *A blinded observational study in 25,505 women. Across the seven fasting glucose septiles, LGA rose 5.3% → 26.3% and cord C-peptide >90th percentile 3.7% → 32.4%. Established that this is a continuous gradient rather than a threshold.*
   `[→ the primary standard for calibrating P(LGA), scenario S9]`

2. Crowther CA, et al. **Effect of treatment of gestational diabetes mellitus on pregnancy outcomes.** N Engl J Med 2005;352:2477-86. (ACHOIS)
   <https://pubmed.ncbi.nlm.nih.gov/15951574/>
   *Treatment reduced the serious perinatal composite outcome from 4% to 1%.* `[→ treatment effect size]`

3. Landon MB, et al. **A multicenter, randomized trial of treatment for mild gestational diabetes.** N Engl J Med 2009;361:1339-48. (MFMU)
   <https://pubmed.ncbi.nlm.nih.gov/19797280/>
   *In mild GDM, LGA 14.5→7.1%, macrosomia 14.3→5.9%, shoulder dystocia 4.0→1.5%, pre-eclampsia 13.6→8.6%, caesarean section 33.8→26.9%.*
   `[→ PSD0/KSD, PCS0/KCS, PPE0/KPEG/KPEB were solved analytically and directly from these five pairs]`

4. Langer O, et al. **A comparison of glyburide and insulin in women with gestational diabetes mellitus.** N Engl J Med 2000;343:1134-8.
   <https://pubmed.ncbi.nlm.nih.gov/11036118/>
   *The original study that made glyburide widely used as the first alternative. Later studies weakened this conclusion (see §5).*

5. Rowan JA, et al. **Metformin versus insulin for the treatment of gestational diabetes.** N Engl J Med 2008;358:2003-15. (MiG)
   <https://pubmed.ncbi.nlm.nih.gov/18463376/>
   *The composite neonatal outcome was equivalent. The key figure: **46% of the metformin arm needed supplemental insulin**.*
   `[→ scenario S8 is an attempt to reproduce this 46% as a prediction rather than as an input]`

6. Sénat MV, et al. **Effect of Glyburide vs Subcutaneous Insulin on Perinatal Complications Among Women With Gestational Diabetes: A Randomized Clinical Trial.** JAMA 2018;319:1773-1780.
   <https://pubmed.ncbi.nlm.nih.gov/29715355/>
   *Failed to demonstrate non-inferiority. Neonatal hypoglycaemia was the main component of the excess risk.* `[→ the glyburide fetal SUR1 term]`

7. Balsells M, et al. **Glibenclamide, metformin, and insulin for the treatment of gestational diabetes: a systematic review and meta-analysis.** BMJ 2015;350:h102.
   <https://pubmed.ncbi.nlm.nih.gov/25609400/>
   *Glyburide against insulin: birth weight +109 g, macrosomia RR 2.62, neonatal hypoglycaemia RR 2.04.*

8. Feig DS, et al. **Metformin in women with type 2 diabetes in pregnancy (MiTy): a multicentre, international, randomised, placebo-controlled trial.** Lancet Diabetes Endocrinol 2020;8:834-844.
   <https://pubmed.ncbi.nlm.nih.gov/32946820/>
   *Adding metformin to insulin reduced insulin requirement, weight gain, and macrosomia, but with a signal of increased SGA.*

9. Feig DS, et al. **Continuous glucose monitoring in pregnant women with type 1 diabetes (CONCEPTT): a multicentre international randomised controlled trial.** Lancet 2017;390:2347-2359.
   <https://pubmed.ncbi.nlm.nih.gov/28923465/>
   `[→ the TIR target and the definitions of the CGM indices]`

10. Hillier TA, et al. **A Pragmatic, Randomized Clinical Trial of Gestational Diabetes Screening.** N Engl J Med 2021;384:895-904. (ScreenR2GDM)
    <https://pubmed.ncbi.nlm.nih.gov/33704936/>
    *The one-step (IADPSG) approach doubled the diagnosis rate but did not improve perinatal outcomes — the reason the diagnostic cut-point problem has to be modelled.*

11. Stewart ZA, et al. **Closed-Loop Insulin Delivery during Pregnancy in Women with Type 1 Diabetes.** N Engl J Med 2016;375:644-54.
    <https://pubmed.ncbi.nlm.nih.gov/27532830/>

12. Donovan LE, et al. **Closed-Loop Insulin Delivery in Type 1 Diabetes in Pregnancy: The CIRCUIT Randomized Clinical Trial.** JAMA 2025;334:2176-2185.
    <https://pubmed.ncbi.nlm.nih.gov/41134589/>

---

## 2. Mechanisms of gestational insulin resistance

The basis for the model's `CID` (counter-insulin index) → `SIREL` block.

13. Barbour LA, et al. **Cellular mechanisms for insulin resistance in normal pregnancy and gestational diabetes.** Diabetes Care 2007;30 Suppl 2:S112-9.
    <https://pubmed.ncbi.nlm.nih.gov/17596458/>
    *A comprehensive review arguing that post-receptor defects (reduced IRS-1 protein, increased p85, reduced GLUT4 translocation) are the core of it.*
    `[→ the basis for modelling SIREL as a post-receptor gain]`

14. Kirwan JP, et al. **TNF-alpha is a predictor of insulin resistance in human pregnancy.** Diabetes 2002;51:2207-13.
    <https://pubmed.ncbi.nlm.nih.gov/12086951/>
    *TNF-α correlates inversely with insulin sensitivity more strongly than hPL, cortisol, oestradiol, prolactin, or leptin.*
    `[→ why STNF carries the largest weight per unit change]`

15. Catalano PM, et al. **Longitudinal changes in glucose metabolism during pregnancy in obese women with normal glucose tolerance and gestational diabetes mellitus.** Am J Obstet Gynecol 1999;180:903-16.
    <https://pubmed.ncbi.nlm.nih.gov/10203659/>
    `[→ the target of a 50-60% fall in SI during pregnancy]`

16. Catalano PM, et al. **Carbohydrate metabolism during pregnancy in control subjects and women with gestational diabetes.** Am J Physiol 1993;264:E60-7.
    <https://pubmed.ncbi.nlm.nih.gov/8430789/>
    `[→ changes in EGP and in basal glucose turnover during pregnancy]`

17. Catalano PM, et al. **Longitudinal changes in insulin release and insulin resistance in nonobese pregnant women.** Am J Obstet Gynecol 1991;165:1667-72.
    <https://pubmed.ncbi.nlm.nih.gov/1750458/>

18. Ryan EA, et al. **Insulin action during pregnancy. Studies with the euglycemic clamp technique.** Diabetes 1985;34:380-9.
    <https://pubmed.ncbi.nlm.nih.gov/3882502/>
    *The classic clamp study of insulin resistance in pregnancy.*

19. Barbour LA, et al. **Human placental growth hormone causes severe insulin resistance in transgenic mice.** Am J Obstet Gynecol 2002;186:512-7.
    <https://pubmed.ncbi.nlm.nih.gov/11904616/>
    `[→ the basis for including GH-V alongside hPL in the placental drive]`

20. Barbour LA, et al. **Human placental growth hormone increases expression of the p85 regulatory unit of phosphatidylinositol 3-kinase and triggers severe insulin resistance in skeletal muscle.** Endocrinology 2004;145:1144-50.
    <https://pubmed.ncbi.nlm.nih.gov/14633976/>

21. Masuyama H, Hiramatsu Y. **Potential role of estradiol and progesterone in insulin resistance through constitutive androstane receptor.** J Mol Endocrinol 2011;47:229-39.
    <https://pubmed.ncbi.nlm.nih.gov/21768169/>
    `[→ the SPROG weighting]`

22. González C, et al. **Regulation of insulin receptor substrate-1 in the liver, skeletal muscle and adipose tissue of rats throughout pregnancy.** Gynecol Endocrinol 2003;17:187-97.
    <https://pubmed.ncbi.nlm.nih.gov/12857426/>
    `[→ separation of hepatic from peripheral resistance (HFR = 0.8)]`

---

## 3. β-cell adaptation and its failure

**This section is the basis for the central hypothesis of the model (BCAP).**

23. Butler AE, et al. **Adaptive changes in pancreatic beta cell fractional area and beta cell turnover in human pregnancy.** Diabetologia 2010;53:2167-76.
    <https://pubmed.ncbi.nlm.nih.gov/20523966/>
    *In humans the increase in β-cell mass is only about 1.4-fold — unlike the 2-3-fold of rodents. Human adaptation is centred on **function rather than mass**.*
    `[→ why the functional component was separated out through a BCM ceiling of 2.6 and KSENS (a leftward shift of the curve)]`

24. Kim H, et al. **Serotonin regulates pancreatic beta cell mass during pregnancy.** Nat Med 2010;16:804-8.
    <https://pubmed.ncbi.nlm.nih.gov/20581837/>
    *Lactogen → Tph1 → serotonin → HTR2B (proliferation)/HTR3A (a lowered secretory threshold).*
    `[→ the molecular basis for the ADAPT · KSENS pathway]`

25. Schraenen A, et al. **Placental lactogens induce serotonin biosynthesis in a subset of mouse beta cells during pregnancy.** Diabetologia 2010;53:2589-99.
    <https://pubmed.ncbi.nlm.nih.gov/20938637/>

26. Vasavada RC, et al. **Targeted expression of placental lactogen in the beta cells of transgenic mice results in beta cell proliferation, islet mass augmentation, and hypoglycemia.** J Biol Chem 2000;275:15399-406.
    <https://pubmed.ncbi.nlm.nih.gov/10809775/>
    `[→ the LACTDR = HPL/(HPL+KMH) driving term]`

27. Billestrup N, Nielsen JH. **The stimulatory effect of growth hormone, prolactin, and placental lactogen on beta-cell proliferation is not mediated by insulin-like growth factor-I.** Endocrinology 1991;129:883-8.
    <https://pubmed.ncbi.nlm.nih.gov/1677331/>

28. Qiao L, et al. **Adiponectin Promotes Maternal β-Cell Expansion Through Placental Lactogen Expression.** Diabetes 2021;70:132-142.
    <https://pubmed.ncbi.nlm.nih.gov/33087456/>
    *Adiponectin → placental lactogen → β-cell expansion. Low adiponectin may be causal rather than simply a marker.*
    `[→ the reason ADIPO is held as a state variable; the current model connects only the SI pathway and leaves this β-cell pathway unimplemented (stated as a limitation)]`

29. Buchanan TA, Xiang AH. **Gestational diabetes mellitus.** J Clin Invest 2005;115:485-91.
    <https://pubmed.ncbi.nlm.nih.gov/15765129/>
    *The review that framed GDM as "a chronic β-cell defect revealed by the stress test that is pregnancy".*
    `[→ the entire central thesis of the model]`

30. Powe CE, et al. **Heterogeneous Contribution of Insulin Sensitivity and Secretion Defects to Gestational Diabetes Mellitus.** Diabetes Care 2016;39:1052-5.
    <https://pubmed.ncbi.nlm.nih.gov/27208340/>
    *Decomposes GDM into an insulin-resistance-predominant, a secretion-deficient, and a mixed form. Phenotype and treatment response differ by subtype.*
    `[→ in this model the subtype is not assigned but emerges from BCAP×BMI — the subtype coordinates on tab 1]`

31. Xiang AH, et al. **Longitudinal changes in insulin sensitivity and beta cell function between women with and without a history of gestational diabetes mellitus.** Diabetologia 2013;56:2753-60.
    <https://pubmed.ncbi.nlm.nih.gov/24030069/>
    `[→ the postpartum decline in β-cells (KDECL) and the disposition index trajectory]`

---

## 4. Placental transport and the fetal compartment

32. Freinkel N. **Banting Lecture 1980. Of pregnancy and progeny.** Diabetes 1980;29:1023-35.
    <https://pubmed.ncbi.nlm.nih.gov/7002669/>
    *"Fuel-mediated teratogenesis" — the proposition that amino acids and lipids, not glucose alone, are fuels for fetal programming.*
    `[→ the basis for the BFAT (maternal FFA → fetal fat) term]`

33. Freinkel N, et al. **Diabetic embryopathy and fuel-mediated organ teratogenesis: lessons from animal models.** Horm Metab Res 1988;20:463-75.
    <https://pubmed.ncbi.nlm.nih.gov/3053387/>

34. Hahn T, et al. **Sustained hyperglycemia in vitro down-regulates the GLUT1 glucose transport system of cultured human term placental trophoblast: a mechanism to protect fetal development?** FASEB J 1998;12:1221-31.
    <https://pubmed.ncbi.nlm.nih.gov/9737725/>
    *A protective mechanism by which hyperglycaemia downregulates GLUT1 — the current model **does not implement** this adaptation, and so may overestimate fetal exposure at extreme hyperglycaemia.*
    `[→ stated as a limitation]`

35. Takata K, et al. **Immunolocalization of glucose transporter GLUT1 in the rat placental barrier.** Cell Tissue Res 1994;276:411-8.
    <https://pubmed.ncbi.nlm.nih.gov/8062336/>

36. Barta E, Drugan A. **A theoretical model of glucose transport suggests symmetric GLUT1 characteristics at placental membranes.** J Membr Biol 2014;247:685-94.
    <https://pubmed.ncbi.nlm.nih.gov/24894722/>
    `[→ the basis for modelling KTRF as symmetric bidirectional diffusion]`

37. Catalano PM, et al. **Increased fetal adiposity: a very sensitive marker of abnormal in utero development.** Am J Obstet Gynecol 2003;189:1698-704.
    <https://pubmed.ncbi.nlm.nih.gov/14710101/>
    *At the same birth weight, GDM neonates carry more body fat — direct evidence that the overgrowth is asymmetric.*
    `[→ the key basis for AFAT (0.55) ≫ AIGF (0.20); normal term fat 13.3%]`

38. Catalano PM, et al. **Fetuses of obese mothers develop insulin resistance in utero.** Diabetes Care 2009;32:1076-80.
    <https://pubmed.ncbi.nlm.nih.gov/19460915/>

39. Schaefer-Graf UM, et al. **Differences in the implications of maternal lipids on fetal metabolism and growth between gestational diabetes mellitus and control pregnancies.** Diabet Med 2011;28:1053-9.
    <https://pubmed.ncbi.nlm.nih.gov/21658120/>
    `[→ the BFAT coefficient]`

40. Kemball ML, et al. **Neonatal hypoglycaemia in infants of diabetic mothers given sulphonylurea drugs in pregnancy.** Arch Dis Child 1970;45:696-701.
    <https://pubmed.ncbi.nlm.nih.gov/5477685/>
    *An observation from 50 years ago that fetal exposure to sulphonylureas causes neonatal hypoglycaemia — a foreshadowing of Sénat's 2018 result.*

41. Haworth JC, et al. **Prognosis of infants of diabetic mothers in relation to neonatal hypoglycaemia.** Dev Med Child Neurol 1976;18:471-9.
    <https://pubmed.ncbi.nlm.nih.gov/955311/>

42. Agrawal RK, et al. **Neonatal hypoglycaemia in infants of diabetic mothers.** J Paediatr Child Health 2000;36:354-6.
    <https://pubmed.ncbi.nlm.nih.gov/10940170/>
    `[→ the PNH0 / PNHMAX range]`

---

## 5. Drug PK/PD and placental passage

**The only axis on which the three drugs genuinely differ.**

### Metformin (crosses freely)

43. Eyal S, et al. **Pharmacokinetics of metformin during pregnancy.** Drug Metab Dispos 2010;38:833-40.
    <https://pubmed.ncbi.nlm.nih.gov/20118196/>
    `[→ CLM · CLMP (renal clearance +25% in pregnancy)]`

44. Charles B, et al. **Population pharmacokinetics of metformin in late pregnancy.** Ther Drug Monit 2006;28:67-72.
    <https://pubmed.ncbi.nlm.nih.gov/16418696/>
    `[→ the VMC/VMP/QM two-compartment structure]`

45. Abduljalil K, et al. **Prediction of Maternal and Fetal Acyclovir, Emtricitabine, Lamivudine, and Metformin Concentrations during Pregnancy Using a Physiologically Based Pharmacokinetic Modeling Approach.** Clin Pharmacokinet 2022;61:725-748.
    <https://pubmed.ncbi.nlm.nih.gov/35067869/>
    *A maternal-fetal PBPK. What KPLM/CLMF in this QSP model has to reproduce.*

46. Tiley JB, et al. **Comparison of Metformin PBPK Models Incorporating Placental Transfer to Predict Fetal and Maternal Exposure.** CPT Pharmacometrics Syst Pharmacol 2026;15:e70136.
    <https://pubmed.ncbi.nlm.nih.gov/41289433/>
    `[→ the umbilical:maternal ≈ 1.0 target; the model gives KPLM/(KPLM+CLMF) = 0.94]`

47. Gu X, et al. **Transplacental transfer of metformin and interaction of metformin with the uptake transporters of placental trophoblast cells.** Eur J Pharm Sci 2025;212:107161.
    <https://pubmed.ncbi.nlm.nih.gov/40494429/>
    `[→ OCT3/PMAT-mediated uptake]`

48. Sheng B, et al. **Short-term neonatal outcomes in women with gestational diabetes treated using metformin versus insulin: a systematic review and meta-analysis of randomized controlled trials.** Acta Diabetol 2023;60:595-608.
    <https://pubmed.ncbi.nlm.nih.gov/36593391/>

### Glyburide (crosses; BCRP-effluxed)

49. Hebert MF, et al. **Are we optimizing gestational diabetes treatment with glyburide? The pharmacologic basis for better clinical practice.** Clin Pharmacol Ther 2009;85:607-14.
    <https://pubmed.ncbi.nlm.nih.gov/19295505/>
    *Glyburide clearance roughly doubles in pregnancy — underdosing and fetal exposure are a problem at the same time.*
    `[→ CLGP = 1.0 (2-fold)]`

50. Gedeon C, et al. **Breast cancer resistance protein: mediating the trans-placental transfer of glyburide across the human placenta.** Placenta 2008;29:39-43.
    <https://pubmed.ncbi.nlm.nih.gov/17923155/>
    `[→ defines KPLG as the net transport with BCRP efflux subtracted]`

51. Kraemer J, et al. **Perfusion studies of glyburide transfer across the human placenta: implications for fetal safety.** Am J Obstet Gynecol 2006;195:270-4.
    <https://pubmed.ncbi.nlm.nih.gov/16579925/>

52. Elliott BD, et al. **Insignificant transfer of glyburide occurs across the human placenta.** Am J Obstet Gynecol 1991;165:807-12.
    <https://pubmed.ncbi.nlm.nih.gov/1951536/>
    *The original paper that reached the opposite conclusion — later perfusion and cord blood studies overturned it. Where the literature conflicts, it is honest to state which side the model took: this model accepts passage and sets cord:maternal ≈ 0.7.*

53. Yu DQ, et al. **Glycemic control and neonatal outcomes in women with gestational diabetes mellitus treated using glyburide, metformin, or insulin: a pairwise and network meta-analysis.** BMC Endocr Disord 2021;21:199.
    <https://pubmed.ncbi.nlm.nih.gov/34641848/>

### Insulin (does not cross)

54. Athanasiadou KI, et al. **Safety and efficacy of insulin detemir versus NPH in the treatment of diabetes during pregnancy: Systematic review and meta-analysis of randomized controlled trials.** Diabetes Res Clin Pract 2022;190:110020.
    <https://pubmed.ncbi.nlm.nih.gov/35878788/>
    `[→ KAB (basal analogue absorption)]`

55. Koren R, Toledano Y, Hod M. **The use of insulin detemir during pregnancy: a safety evaluation.** Expert Opin Drug Saf 2015;14:593-9.
    <https://pubmed.ncbi.nlm.nih.gov/25731934/>

---

## 6. Long-term maternal outcomes

56. Bellamy L, et al. **Type 2 diabetes mellitus after gestational diabetes: a systematic review and meta-analysis.** Lancet 2009;373:1773-9.
    <https://pubmed.ncbi.nlm.nih.gov/19465232/>
    *A relative risk of 7.43.* `[→ KDI = 3.8 was back-calculated from this RR; at DI/DIREF = 0.45, RR ≈ 8]`

57. Vounzoulaki E, et al. **Progression to type 2 diabetes in women with a known history of gestational diabetes: systematic review and meta-analysis.** BMJ 2020;369:m1361.
    <https://pubmed.ncbi.nlm.nih.gov/32404325/>
    `[→ the five-year cumulative incidence]`

58. Gunderson EP, et al. **Lactation and Progression to Type 2 Diabetes Mellitus After Gestational Diabetes Mellitus: A Prospective Cohort Study.** Ann Intern Med 2015;163:889-98.
    <https://pubmed.ncbi.nlm.nih.gov/26595611/>
    `[→ LACT / KLACT = 0.45]`

59. Schwartz N, et al. **The prevalence of gestational diabetes mellitus recurrence--effect of ethnicity and parity: a metaanalysis.** Am J Obstet Gynecol 2015;213:310-7.
    <https://pubmed.ncbi.nlm.nih.gov/25757637/>

60. Pan Y, et al. **Gestational diabetes mellitus recurrence rate and risk factors: a systematic review and meta-analysis.** Diabetes Res Clin Pract 2025;230:112949.
    <https://pubmed.ncbi.nlm.nih.gov/41130422/>

61. Damm P, et al. **Gestational diabetes mellitus and long-term consequences for mother and offspring: a view from Denmark.** Diabetologia 2016;59:1396-1399.
    <https://pubmed.ncbi.nlm.nih.gov/27174368/>

62. Nouhjah S, et al. **Postpartum screening practices, progression to abnormal glucose tolerance and its related risk factors in Asian women with a known history of gestational diabetes: A systematic review and meta-analysis.** Diabetes Metab Syndr 2017;11 Suppl 2:S703-S712.
    <https://pubmed.ncbi.nlm.nih.gov/28571777/>

---

## 7. Offspring outcomes and developmental programming

63. Lowe WL Jr, et al. **Association of Gestational Diabetes With Maternal Disorders of Glucose Metabolism and Childhood Adiposity.** JAMA 2018;320:1005-1016. (HAPO Follow-Up Study)
    <https://pubmed.ncbi.nlm.nih.gov/30208453/>

64. Hillier TA, et al. **Childhood obesity and metabolic imprinting: the ongoing effects of maternal hyperglycemia.** Diabetes Care 2007;30:2287-92.
    <https://pubmed.ncbi.nlm.nih.gov/17519427/>

65. Dabelea D, Crume T. **Maternal environment and the transgenerational cycle of obesity and diabetes.** Diabetes 2011;60:1849-55.
    <https://pubmed.ncbi.nlm.nih.gov/21709280/>

66. Rowan JA, et al. **Metformin in gestational diabetes: the offspring follow-up (MiG TOFU): body composition at 2 years of age.** Diabetes Care 2011;34:2279-84.
    <https://pubmed.ncbi.nlm.nih.gov/21949222/>
    *Metformin-exposed children had more subcutaneous fat at age two.*
    `[→ the model does not predict this result; the fetal metformin term of dxdt_FATF is deliberately left at 0 and inactive, to expose that gap]`

67. Rowan JA, et al. **Metformin in Gestational Diabetes The Offspring Follow Up (MiGTOFU): Associations between maternal characteristics and size and adiposity of boys and girls at nine years.** Aust N Z J Obstet Gynaecol 2023;63:825-828.
    <https://pubmed.ncbi.nlm.nih.gov/37469163/>

68. Plagemann A, Harder T. **Fuel-mediated teratogenesis and breastfeeding.** Diabetes Care 2011;34:779-81.
    <https://pubmed.ncbi.nlm.nih.gov/21357365/>

---

## 8. Diagnosis, screening, guidelines

69. International Association of Diabetes and Pregnancy Study Groups Consensus Panel. **International association of diabetes and pregnancy study groups recommendations on the diagnosis and classification of hyperglycemia in pregnancy.** Diabetes Care 2010;33:676-82.
    <https://pubmed.ncbi.nlm.nih.gov/20190296/>
    *The IADPSG criteria: fasting 5.1 / 1 h 10.0 / 2 h 8.5 mmol/L (92 / 180 / 153 mg/dL).*

70. American Diabetes Association. **13. Management of Diabetes in Pregnancy: Standards of Medical Care in Diabetes-2018.** Diabetes Care 2018;41:S137-S143.
    <https://pubmed.ncbi.nlm.nih.gov/29222384/>
    `[→ the treatment targets of fasting <95 and 1 h <140 mg/dL]`

71. **ACOG Practice Bulletin No. 190: Gestational Diabetes Mellitus.** Obstet Gynecol 2018;131:e49-e64.
    <https://pubmed.ncbi.nlm.nih.gov/29370047/>

72. Brady M, et al. **One-Step Compared With Two-Step Gestational Diabetes Screening and Pregnancy Outcomes: A Systematic Review and Meta-analysis.** Obstet Gynecol 2022;140:712-723.
    <https://pubmed.ncbi.nlm.nih.gov/36201772/>

---

## 9. Biomarkers

73. Ghosh P, et al. **Plasma Glycated CD59, a Novel Biomarker for Detection of Pregnancy-Induced Glucose Intolerance.** Diabetes Care 2017;40:981-984.
    <https://pubmed.ncbi.nlm.nih.gov/28450368/>

74. Williams MA, et al. **Plasma adiponectin concentrations in early pregnancy and subsequent risk of gestational diabetes mellitus.** J Clin Endocrinol Metab 2004;89:2306-11.
    <https://pubmed.ncbi.nlm.nih.gov/15126557/>
    `[→ the basis for holding ADIPO as a first-trimester predictor]`

75. Ye Y, et al. **Adiponectin, leptin, and leptin/adiponectin ratio with risk of gestational diabetes mellitus: A prospective nested case-control study among Chinese women.** Diabetes Res Clin Pract 2022;191:110039.
    <https://pubmed.ncbi.nlm.nih.gov/35985429/>

76. Durnwald C, et al. **Continuous Glucose Monitoring Profiles in Pregnancies With and Without Gestational Diabetes Mellitus.** Diabetes Care 2024;47:1333-1341.
    <https://pubmed.ncbi.nlm.nih.gov/38701400/>
    `[→ the reference values for TIR and mean glucose on the CGM tab]`

77. García-Moreno RM, et al. **Efficacy of continuous glucose monitoring on maternal and neonatal outcomes in gestational diabetes mellitus: a systematic review and meta-analysis of randomized clinical trials.** Diabet Med 2022;39:e14703.
    <https://pubmed.ncbi.nlm.nih.gov/34564868/>

78. Benhalima K, et al. **Application of continuous glucose monitoring and automated insulin delivery technologies for pregnant women with type 1, type 2, or gestational diabetes: an international consensus statement.** Lancet Diabetes Endocrinol 2026;14:157-177.
    <https://pubmed.ncbi.nlm.nih.gov/41421368/>

---

## 10. Prevention: lifestyle, exercise, supplements

79. Davenport MH, et al. **Prenatal exercise for the prevention of gestational diabetes mellitus and hypertensive disorders of pregnancy: a systematic review and meta-analysis.** Br J Sports Med 2018;52:1367-1375.
    <https://pubmed.ncbi.nlm.nih.gov/30337463/>
    `[→ EXEFF = 0.20 (150 minutes a week or more)]`

80. Tsironikos GI, et al. **Effectiveness of exercise intervention during pregnancy on high-risk women for gestational diabetes mellitus prevention: A meta-analysis of published RCTs.** PLoS One 2022;17:e0272711.
    <https://pubmed.ncbi.nlm.nih.gov/35930592/>

81. Mashayekh-Amiri S, et al. **Myo-inositol supplementation for prevention of gestational diabetes mellitus in overweight and obese pregnant women: a systematic review and meta-analysis.** Diabetol Metab Syndr 2022;14:93.
    <https://pubmed.ncbi.nlm.nih.gov/35794663/>

---

## 11. Genetics of GDM susceptibility

82. Shan D, et al. **MTNR1B rs1387153 Polymorphism and Risk of Gestational Diabetes Mellitus: Meta-Analysis and Trial Sequential Analysis.** Public Health Genomics 2023;26:201-211.
    <https://pubmed.ncbi.nlm.nih.gov/37980891/>
    `[→ the basis for holding BCAP as a patient-level parameter that includes a genetic component]`

83. Zhang Y, et al. **MTNR1B gene variations and high pre-pregnancy BMI increase gestational diabetes mellitus risk in Chinese women.** Gene 2024;894:148023.
    <https://pubmed.ncbi.nlm.nih.gov/38007162/>
    *The clinical counterpart of the BCAP (genetic) × BMI (resistance) interaction — the subtype coordinates on tab 1 of the model are exactly this two-dimensional plane.*

84. Huang LT, et al. **Adiponectin gene polymorphisms and risk of gestational diabetes mellitus: A meta-analysis.** World J Clin Cases 2019;7:572-584.
    <https://pubmed.ncbi.nlm.nih.gov/30863757/>

---

## 12. Shared biology with preeclampsia

85. Elgazzaz M, Lazartigues E. **Implications of pregnancy on cardiometabolic disease risk: preeclampsia and gestational diabetes.** Am J Physiol Cell Physiol 2024;327:C646-C660.
    <https://pubmed.ncbi.nlm.nih.gov/39010840/>
    *GDM and pre-eclampsia share their upstream biology (endothelial dysfunction, inflammation). sFlt-1/PlGF is drawn on the map of this model but **is not implemented as a differential equation**, and pre-eclampsia is handled by a regression equation only — stated as a limitation.*

86. Aziz F, et al. **Gestational diabetes mellitus, hypertension, and dyslipidemia as the risk factors of preeclampsia.** Sci Rep 2024;14:6182.
    <https://pubmed.ncbi.nlm.nih.gov/38486097/>

87. Verlohren S, Dröge LA. **Clinical interpretation and implementation of the sFlt-1/PlGF ratio in the prediction, diagnosis and management of preeclampsia.** Pregnancy Hypertens 2022;27:42-50.
    <https://pubmed.ncbi.nlm.nih.gov/34915395/>

---

## Where the literature disagrees

The model **chose** one side at the following three points. It is written down here because, if the choice is not stated,
the parameter looks like a fact.

| Issue | Conflicting evidence | This model's choice |
|---|---|---|
| Placental passage of glyburide | Elliott 1991 (#52), "negligible passage", vs the perfusion/BCRP studies of Kraemer 2006 (#51) and Gedeon 2008 (#50) | passage accepted, cord:maternal ≈ 0.70 |
| One-step vs two-step screening | the IADPSG (#69) recommendation vs Hillier NEJM 2021 (#10), no improvement in perinatal outcome | both criteria are output as observed values, and neither is treated as the right answer |
| Long-term effects of metformin on the offspring | short-term equivalence in MiG (#5) vs increased fat in MiG-TOFU (#66, #67) | fetal exposure is computed but the growth effect is **fixed at 0** — putting in an effect with no mechanism behind it would be a fudge |

---

## Caveat

The purpose of this reference list is to make the **provenance** of the model parameters traceable.
Some individual parameters are values reported as such in the paper concerned, and others were
back-calculated by hand from reported steady-state observations; the two are mixed together. Which is which is
recorded item by item in the CALIBRATION section of `gdm_mrgsolve_model_en.R`.
The model has never been fitted to individual patient data.
