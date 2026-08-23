# Kidney Transplant Rejection — References

`ktx_qsp_model.dot` (mechanistic map), `ktx_mrgsolve_model.R` (ODE model),
`ktx_shiny_app.R` (dashboard) — the literature used for their structure, parameters and validation targets.

**Every PMID was actually looked up through the NCBI E-utilities and confirmed, and the title,
journal and year were taken verbatim from the PubMed record** (these are not citations from
memory). Every link takes the form `https://pubmed.ncbi.nlm.nih.gov/<PMID>/`.

---

## A. Epidemiology · Long-Term Survival · Causes of Graft Loss

> Where long-term graft survival currently stands, and the actual distribution of causes of graft loss. What the model's final endpoints (death-censored graft survival · eGFR slope) aim at.

1. Hariharan S et al. **Long-Term Survival after Kidney Transplantation.** *N Engl J Med* 2021. [PMID 34407344](https://pubmed.ncbi.nlm.nih.gov/34407344/)
2. Coemans M et al. **Analyses of the short- and long-term graft survival after kidney transplantation in Europe between 1986 and 2015.** *Kidney Int* 2018. [PMID 30049474](https://pubmed.ncbi.nlm.nih.gov/30049474/)
3. Lamb KE et al. **Long-term renal allograft survival in the United States: a critical reappraisal.** *Am J Transplant* 2011. [PMID 20973913](https://pubmed.ncbi.nlm.nih.gov/20973913/)
4. Sellarés J et al. **Understanding the causes of kidney transplant failure: the dominant role of antibody-mediated rejection and nonadherence.** *Am J Transplant* 2012. [PMID 22081892](https://pubmed.ncbi.nlm.nih.gov/22081892/)
5. Gaston RS et al. **Evidence for antibody-mediated injury as a major determinant of late kidney allograft failure.** *Transplantation* 2010. [PMID 20463643](https://pubmed.ncbi.nlm.nih.gov/20463643/)
6. El-Zoghby ZM et al. **Identifying specific causes of kidney allograft loss.** *Am J Transplant* 2009. [PMID 19191769](https://pubmed.ncbi.nlm.nih.gov/19191769/)
7. Mayrdorfer M et al. **Exploring the Complexity of Death-Censored Kidney Allograft Failure.** *J Am Soc Nephrol* 2021. [PMID 33883251](https://pubmed.ncbi.nlm.nih.gov/33883251/)
8. Wekerle T et al. **Strategies for long-term preservation of kidney graft function.** *Lancet* 2017. [PMID 28561006](https://pubmed.ncbi.nlm.nih.gov/28561006/)
9. Loupy A et al. **Prediction system for risk of allograft loss in patients receiving kidney transplants: international derivation and validation study.** *BMJ* 2019. [PMID 31530561](https://pubmed.ncbi.nlm.nih.gov/31530561/)
10. Clayton PA et al. **Relationship between eGFR Decline and Hard Outcomes after Kidney Transplants.** *J Am Soc Nephrol* 2016. [PMID 27059513](https://pubmed.ncbi.nlm.nih.gov/27059513/)

## B. Alloantigen Recognition · HLA Mismatch · Eplets

> How alloantigen is recognised — the direct, indirect and semi-direct pathways, and HLA mismatch and eplet load. The basis of the model's AG (antigen presentation load) and BEPI (B-cell epitope load) parameters.

11. Lechler RI et al. **Restoration of immunogenicity to passenger cell-depleted kidney allografts by the addition of donor strain dendritic cells.** *J Exp Med* 1982. [PMID 7033437](https://pubmed.ncbi.nlm.nih.gov/7033437/)
12. Ali JM et al. **Allorecognition pathways in transplant rejection and tolerance.** *Transplantation* 2013. [PMID 23715047](https://pubmed.ncbi.nlm.nih.gov/23715047/)
13. Harper SJ et al. **CD8 T-cell recognition of acquired alloantigen promotes acute allograft rejection.** *Proc Natl Acad Sci U S A* 2015. [PMID 26420874](https://pubmed.ncbi.nlm.nih.gov/26420874/)
14. Opelz G et al. **Effect of human leukocyte antigen compatibility on kidney graft survival: comparative analysis of two decades.** *Transplantation* 2007. [PMID 17667803](https://pubmed.ncbi.nlm.nih.gov/17667803/)
15. Duquesnoy RJ et al. **HLAmatchmaker: a molecularly based algorithm for histocompatibility determination. IV. An alternative strategy to increase the number of compatible donors for highly sensitized patients.** *Transplantation* 2003. [PMID 12660520](https://pubmed.ncbi.nlm.nih.gov/12660520/)
16. Wiebe C et al. **Class II HLA epitope matching-A strategy to minimize de novo donor-specific antibody development and improve outcomes.** *Am J Transplant* 2013. [PMID 24164958](https://pubmed.ncbi.nlm.nih.gov/24164958/)
17. Senev A et al. **Eplet Mismatch Load and De Novo Occurrence of Donor-Specific Anti-HLA Antibodies, Rejection, and Graft Failure after Kidney Transplantation: An Observational Cohort Study.** *J Am Soc Nephrol* 2020. [PMID 32764139](https://pubmed.ncbi.nlm.nih.gov/32764139/)
18. Wiebe C et al. **HLA-DR/DQ molecular mismatch: A prognostic biomarker for primary alloimmunity.** *Am J Transplant* 2019. [PMID 30414349](https://pubmed.ncbi.nlm.nih.gov/30414349/)

## C. T-Cell-Mediated Rejection (TCMR)

> The mechanism, pathology and treatment response of T-cell-mediated rejection (TCMR). The model's TN→TACT→TEFF→TINF→TUB axis.

19. Halloran PF **T cell-mediated rejection of kidney transplants: a personal viewpoint.** *Am J Transplant* 2010. [PMID 20346061](https://pubmed.ncbi.nlm.nih.gov/20346061/)
20. Nankivell BJ et al. **Rejection of the kidney allograft.** *N Engl J Med* 2010. [PMID 20925547](https://pubmed.ncbi.nlm.nih.gov/20925547/)
21. Bouatou Y et al. **Response to treatment and long-term outcomes in kidney transplant recipients with acute T cell-mediated rejection.** *Am J Transplant* 2019. [PMID 30748089](https://pubmed.ncbi.nlm.nih.gov/30748089/)
22. Nankivell BJ et al. **The clinical and pathological significance of borderline T cell-mediated rejection.** *Am J Transplant* 2019. [PMID 30501008](https://pubmed.ncbi.nlm.nih.gov/30501008/)
23. Halloran PF et al. **Disappearance of T Cell-Mediated Rejection Despite Continued Antibody-Mediated Rejection in Late Kidney Transplant Recipients.** *J Am Soc Nephrol* 2015. [PMID 25377077](https://pubmed.ncbi.nlm.nih.gov/25377077/)
24. Lakkis FG et al. **Origin and biology of the allogeneic response.** *Cold Spring Harb Perspect Med* 2013. [PMID 23906882](https://pubmed.ncbi.nlm.nih.gov/23906882/)
25. Wood KJ et al. **Regulatory immune cells in transplantation.** *Nat Rev Immunol* 2012. [PMID 22627860](https://pubmed.ncbi.nlm.nih.gov/22627860/)
26. Kawai T et al. **HLA-mismatched renal transplantation without maintenance immunosuppression.** *N Engl J Med* 2008. [PMID 18216355](https://pubmed.ncbi.nlm.nih.gov/18216355/)
27. Sawitzki B et al. **Regulatory cell therapy in kidney transplantation (The ONE Study): a harmonised design and analysis of seven non-randomised, single-arm, phase 1/2A trials.** *Lancet* 2020. [PMID 32446407](https://pubmed.ncbi.nlm.nih.gov/32446407/)

## D. B Cells · Donor-Specific Antibody (DSA) · Antibody-Mediated Rejection (ABMR)

> From B cells through to donor-specific antibody (DSA), and antibody-mediated rejection (ABMR). The model's BN→BGC→PB→LLPC→DSA axis and 'why rituximab fails to lower DSA'.

28. Loupy A et al. **Antibody-Mediated Rejection of Solid-Organ Allografts.** *N Engl J Med* 2018. [PMID 30231232](https://pubmed.ncbi.nlm.nih.gov/30231232/)
29. Wiebe C et al. **Evolution and clinical pathologic correlations of de novo donor-specific HLA antibody post kidney transplant.** *Am J Transplant* 2012. [PMID 22429309](https://pubmed.ncbi.nlm.nih.gov/22429309/)
30. Lefaucheur C et al. **Preexisting donor-specific HLA antibodies predict outcome in kidney transplantation.** *J Am Soc Nephrol* 2010. [PMID 20634297](https://pubmed.ncbi.nlm.nih.gov/20634297/)
31. Loupy A et al. **Complement-binding anti-HLA antibodies and kidney-allograft survival.** *N Engl J Med* 2013. [PMID 24066742](https://pubmed.ncbi.nlm.nih.gov/24066742/)
32. Viglietti D et al. **Value of Donor-Specific Anti-HLA Antibody Monitoring and Characterization for Risk Stratification of Kidney Allograft Loss.** *J Am Soc Nephrol* 2017. [PMID 27493255](https://pubmed.ncbi.nlm.nih.gov/27493255/)
33. Everly MJ et al. **Incidence and impact of de novo donor-specific alloantibody in primary renal allografts.** *Transplantation* 2013. [PMID 23380861](https://pubmed.ncbi.nlm.nih.gov/23380861/)
34. Wiebe C et al. **Rates and determinants of progression to graft failure in kidney allograft recipients with de novo donor-specific antibody.** *Am J Transplant* 2015. [PMID 26096305](https://pubmed.ncbi.nlm.nih.gov/26096305/)
35. Schinstock CA et al. **Recommended Treatment for Antibody-mediated Rejection After Kidney Transplantation: The 2019 Expert Consensus From the Transplantion Society Working Group.** *Transplantation* 2020. [PMID 31895348](https://pubmed.ncbi.nlm.nih.gov/31895348/)
36. Sis B et al. **Endothelial gene expression in kidney transplants with alloantibody indicates antibody-mediated damage despite lack of C4d staining.** *Am J Transplant* 2009. [PMID 19681822](https://pubmed.ncbi.nlm.nih.gov/19681822/)
37. Chong AS et al. **Memory B cells in transplantation.** *Transplantation* 2015. [PMID 25525921](https://pubmed.ncbi.nlm.nih.gov/25525921/)
38. Kamburova EG et al. **A single dose of rituximab does not deplete B cells in secondary lymphoid organs but alters phenotype and function.** *Am J Transplant* 2013. [PMID 23570303](https://pubmed.ncbi.nlm.nih.gov/23570303/)
39. Cecka JM **Calculated PRA (CPRA): the new measure of sensitization for transplant candidates.** *Am J Transplant* 2010. [PMID 19958328](https://pubmed.ncbi.nlm.nih.gov/19958328/)

## E. Complement · NK Cells · Microvascular Inflammation

> The two branches by which DSA actually damages tissue — complement (C1q-C4d-MAC) and complement-independent NK ADCC / ‘missing self’. The model's C4D·NKACT·ENDO·MVI.

40. Stegall MD et al. **Terminal complement inhibition decreases antibody-mediated rejection in sensitized renal transplant recipients.** *Am J Transplant* 2011. [PMID 21942930](https://pubmed.ncbi.nlm.nih.gov/21942930/)
41. Sicard A et al. **Detection of C3d-binding donor-specific anti-HLA antibodies at diagnosis of humoral rejection predicts renal graft loss.** *J Am Soc Nephrol* 2015. [PMID 25125383](https://pubmed.ncbi.nlm.nih.gov/25125383/)
42. Koenig A et al. **Missing self triggers NK cell-mediated chronic vascular rejection of solid organ transplants.** *Nat Commun* 2019. [PMID 31767837](https://pubmed.ncbi.nlm.nih.gov/31767837/)
43. Callemeyn J et al. **Missing Self-Induced Microvascular Rejection of Kidney Allografts: A Population-Based Study.** *J Am Soc Nephrol* 2021. [PMID 34301794](https://pubmed.ncbi.nlm.nih.gov/34301794/)
44. Yazdani S et al. **Natural killer cell infiltration is discriminative for antibody-mediated rejection and predicts outcome after kidney transplantation.** *Kidney Int* 2019. [PMID 30396694](https://pubmed.ncbi.nlm.nih.gov/30396694/)
45. Hidalgo LG et al. **NK cell transcripts and NK cells in kidney biopsies from patients with donor-specific antibodies: evidence for NK cell involvement in antibody-mediated rejection.** *Am J Transplant* 2010. [PMID 20659089](https://pubmed.ncbi.nlm.nih.gov/20659089/)
46. Jane-wit D et al. **Complement membrane attack complexes activate noncanonical NF-κB by forming an Akt+ NIK+ signalosome on Rab5+ endosomes.** *Proc Natl Acad Sci U S A* 2015. [PMID 26195760](https://pubmed.ncbi.nlm.nih.gov/26195760/)
47. Thomas KA et al. **The perfect storm: HLA antibodies, complement, FcγRs, and endothelium in transplant rejection.** *Trends Mol Med* 2015. [PMID 25801125](https://pubmed.ncbi.nlm.nih.gov/25801125/)
48. Eskandary F et al. **A Randomized Trial of Bortezomib in Late Antibody-Mediated Kidney Transplant Rejection.** *J Am Soc Nephrol* 2018. [PMID 29242250](https://pubmed.ncbi.nlm.nih.gov/29242250/)

## F. Banff Classification · Histology · Molecular Diagnostics

> The Banff classification and molecular diagnostics. The basis of the t·i·g·ptc·cg·ci scores and the dd-cfDNA the model outputs.

49. Loupy A et al. **The Banff 2019 Kidney Meeting Report (I): Updates on and clarification of criteria for T cell- and antibody-mediated rejection.** *Am J Transplant* 2020. [PMID 32463180](https://pubmed.ncbi.nlm.nih.gov/32463180/)
50. Naesens M et al. **The Banff 2022 Kidney Meeting Report: Reappraisal of microvascular inflammation and the role of biopsy-based transcript diagnostics.** *Am J Transplant* 2024. [PMID 38032300](https://pubmed.ncbi.nlm.nih.gov/38032300/)
51. Roufosse C et al. **A 2018 Reference Guide to the Banff Classification of Renal Allograft Pathology.** *Transplantation* 2018. [PMID 30028786](https://pubmed.ncbi.nlm.nih.gov/30028786/)
52. Haas M et al. **Banff 2013 meeting report: inclusion of c4d-negative antibody-mediated rejection and antibody-associated arterial lesions.** *Am J Transplant* 2014. [PMID 24472190](https://pubmed.ncbi.nlm.nih.gov/24472190/)
53. Solez K et al. **Banff 07 classification of renal allograft pathology: updates and future directions.** *Am J Transplant* 2008. [PMID 18294345](https://pubmed.ncbi.nlm.nih.gov/18294345/)
54. Halloran PF et al. **Microarray diagnosis of antibody-mediated rejection in kidney transplant biopsies: an international prospective study (INTERCOM).** *Am J Transplant* 2013. [PMID 24119109](https://pubmed.ncbi.nlm.nih.gov/24119109/)
55. Reeve J et al. **Assessing rejection-related disease in kidney transplant biopsies based on archetypal analysis of molecular phenotypes.** *JCI Insight* 2017. [PMID 28614805](https://pubmed.ncbi.nlm.nih.gov/28614805/)
56. Bloom RD et al. **Cell-Free DNA and Active Rejection in Kidney Allografts.** *J Am Soc Nephrol* 2017. [PMID 28280140](https://pubmed.ncbi.nlm.nih.gov/28280140/)
57. Halloran PF et al. **The Trifecta Study: Comparing Plasma Levels of Donor-derived Cell-Free DNA with the Molecular Phenotype of Kidney Transplant Biopsies.** *J Am Soc Nephrol* 2022. [PMID 35058354](https://pubmed.ncbi.nlm.nih.gov/35058354/)
58. Friedewald JJ et al. **Development and clinical validity of a novel blood-based molecular biomarker for subclinical acute rejection following kidney transplant.** *Am J Transplant* 2019. [PMID 29985559](https://pubmed.ncbi.nlm.nih.gov/29985559/)

## G. Chronic Graft Injury · IFTA · Nephron Mass

> Chronic graft injury — IFTA and nephron mass (the Brenner hypothesis). The model's LOOP 2 (hyperfiltration → TGF-β → fibrosis → nephron loss).

59. Nankivell BJ et al. **The natural history of chronic allograft nephropathy.** *N Engl J Med* 2003. [PMID 14668458](https://pubmed.ncbi.nlm.nih.gov/14668458/)
60. Naesens M et al. **Chronic histological damage in early indication biopsies is an independent risk factor for late renal allograft failure.** *Am J Transplant* 2013. [PMID 23136888](https://pubmed.ncbi.nlm.nih.gov/23136888/)
61. Stegall MD et al. **Through a glass darkly: seeking clarity in preventing late kidney transplant failure.** *J Am Soc Nephrol* 2015. [PMID 25097209](https://pubmed.ncbi.nlm.nih.gov/25097209/)
62. Brenner BM et al. **Nephron underdosing: a programmed cause of chronic renal allograft failure.** *Am J Kidney Dis* 1993. [PMID 8494022](https://pubmed.ncbi.nlm.nih.gov/8494022/)
63. Hostetter TH et al. **Hyperfiltration in remnant nephrons: a potentially adverse response to renal ablation.** *Am J Physiol* 1981. [PMID 7246778](https://pubmed.ncbi.nlm.nih.gov/7246778/)
64. Denic A et al. **Single-Nephron Glomerular Filtration Rate in Healthy Adults.** *N Engl J Med* 2017. [PMID 28614683](https://pubmed.ncbi.nlm.nih.gov/28614683/)
65. Amer H et al. **Proteinuria after kidney transplantation, relationship to allograft histology and survival.** *Am J Transplant* 2007. [PMID 17941956](https://pubmed.ncbi.nlm.nih.gov/17941956/)
66. Chapman JR et al. **Chronic renal allograft dysfunction.** *J Am Soc Nephrol* 2005. [PMID 16120819](https://pubmed.ncbi.nlm.nih.gov/16120819/)

## H. Ischaemia-Reperfusion Injury · Delayed Graft Function (DGF) · Donor Quality

> Ischaemia-reperfusion injury and delayed graft function (DGF). The model's AKI initial condition (a function of CIT·DCD·KDPI).

67. Siedlecki A et al. **Delayed graft function in the kidney transplant.** *Am J Transplant* 2011. [PMID 21929642](https://pubmed.ncbi.nlm.nih.gov/21929642/)
68. Yarlagadda SG et al. **Association between delayed graft function and allograft and patient survival: a systematic review and meta-analysis.** *Nephrol Dial Transplant* 2009. [PMID 19103734](https://pubmed.ncbi.nlm.nih.gov/19103734/)
69. Debout A et al. **Each additional hour of cold ischemia time significantly increases the risk of graft failure and mortality following renal transplantation.** *Kidney Int* 2015. [PMID 25229341](https://pubmed.ncbi.nlm.nih.gov/25229341/)
70. Tingle SJ et al. **Normothermic and hypothermic machine perfusion preservation versus static cold storage for deceased donor kidney transplantation.** *Cochrane Database Syst Rev* 2024. [PMID 38979743](https://pubmed.ncbi.nlm.nih.gov/38979743/)
71. Zhao H et al. **Ischemia-Reperfusion Injury Reduces Long Term Renal Graft Survival: Mechanism and Beyond.** *EBioMedicine* 2018. [PMID 29398595](https://pubmed.ncbi.nlm.nih.gov/29398595/)
72. Rao PS et al. **A comprehensive risk quantification score for deceased donor kidneys: the kidney donor risk index.** *Transplantation* 2009. [PMID 19623019](https://pubmed.ncbi.nlm.nih.gov/19623019/)
73. Ponticelli C **Ischaemia-reperfusion injury: a major protagonist in kidney transplantation.** *Nephrol Dial Transplant* 2014. [PMID 24335382](https://pubmed.ncbi.nlm.nih.gov/24335382/)
74. Moers C et al. **Machine perfusion or cold storage in deceased-donor kidney transplantation.** *N Engl J Med* 2009. [PMID 19118301](https://pubmed.ncbi.nlm.nih.gov/19118301/)

## I. Tacrolimus PK/PD · CYP3A5 · TDM · Exposure Variability (IPV)

> Tacrolimus PK/PD · CYP3A5 genotype · TDM · intra-patient variability (IPV). The model's closed-loop TDM controller and the 1.9-fold CYP3A5 clearance.

75. Staatz CE et al. **Clinical pharmacokinetics and pharmacodynamics of tacrolimus in solid organ transplantation.** *Clin Pharmacokinet* 2004. [PMID 15244495](https://pubmed.ncbi.nlm.nih.gov/15244495/)
76. Birdwell KA et al. **Clinical Pharmacogenetics Implementation Consortium (CPIC) Guidelines for CYP3A5 Genotype and Tacrolimus Dosing.** *Clin Pharmacol Ther* 2015. [PMID 25801146](https://pubmed.ncbi.nlm.nih.gov/25801146/)
77. Thervet E et al. **Optimization of initial tacrolimus dose using pharmacogenetic testing.** *Clin Pharmacol Ther* 2010. [PMID 20393454](https://pubmed.ncbi.nlm.nih.gov/20393454/)
78. Woillard JB et al. **Population pharmacokinetic model and Bayesian estimator for two tacrolimus formulations--twice daily Prograf and once daily Advagraf.** *Br J Clin Pharmacol* 2011. [PMID 21284698](https://pubmed.ncbi.nlm.nih.gov/21284698/)
79. Shuker N et al. **Intra-patient variability in tacrolimus exposure: causes, consequences for clinical management.** *Transplant Rev (Orlando)* 2015. [PMID 25687818](https://pubmed.ncbi.nlm.nih.gov/25687818/)
80. Borra LC et al. **High within-patient variability in the clearance of tacrolimus is a risk factor for poor long-term outcome after kidney transplantation.** *Nephrol Dial Transplant* 2010. [PMID 20190242](https://pubmed.ncbi.nlm.nih.gov/20190242/)
81. Gonzales HM et al. **A comprehensive review of the impact of tacrolimus intrapatient variability on clinical outcomes in kidney transplantation.** *Am J Transplant* 2020. [PMID 32406604](https://pubmed.ncbi.nlm.nih.gov/32406604/)
82. Bouamar R et al. **Tacrolimus predose concentrations do not predict the risk of acute rejection after renal transplantation: a pooled analysis from three randomized-controlled clinical trials(†).** *Am J Transplant* 2013. [PMID 23480233](https://pubmed.ncbi.nlm.nih.gov/23480233/)
83. Ekberg H et al. **Reduced exposure to calcineurin inhibitors in renal transplantation.** *N Engl J Med* 2007. [PMID 18094377](https://pubmed.ncbi.nlm.nih.gov/18094377/)
84. Brunet M et al. **Therapeutic Drug Monitoring of Tacrolimus-Personalized Therapy: Second Consensus Report.** *Ther Drug Monit* 2019. [PMID 31045868](https://pubmed.ncbi.nlm.nih.gov/31045868/)
85. Nankivell BJ et al. **Calcineurin inhibitor nephrotoxicity: longitudinal assessment by protocol histology.** *Transplantation* 2004. [PMID 15446315](https://pubmed.ncbi.nlm.nih.gov/15446315/)
86. Naesens M et al. **Calcineurin inhibitor nephrotoxicity.** *Clin J Am Soc Nephrol* 2009. [PMID 19218475](https://pubmed.ncbi.nlm.nih.gov/19218475/)
87. Davis S et al. **Lower tacrolimus exposure and time in therapeutic range increase the risk of de novo donor-specific antibodies in the first year of kidney transplantation.** *Am J Transplant* 2018. [PMID 28925597](https://pubmed.ncbi.nlm.nih.gov/28925597/)
88. Gatault P et al. **Reduction of Extended-Release Tacrolimus Dose in Low-Immunological-Risk Kidney Transplant Recipients Increases Risk of Rejection and Appearance of Donor-Specific Antibodies: A Randomized Study.** *Am J Transplant* 2017. [PMID 27862923](https://pubmed.ncbi.nlm.nih.gov/27862923/)

## J. Mycophenolate PK/PD · Enterohepatic Circulation · Myelosuppression

> Mycophenolate PK/PD and enterohepatic circulation. The model's MPAG-EHC loop and IMPDH inhibition.

89. Staatz CE et al. **Clinical pharmacokinetics and pharmacodynamics of mycophenolate in solid organ transplant recipients.** *Clin Pharmacokinet* 2007. [PMID 17201457](https://pubmed.ncbi.nlm.nih.gov/17201457/)
90. Allison AC **Mechanisms of action of mycophenolate mofetil.** *Lupus* 2005. [PMID 15803924](https://pubmed.ncbi.nlm.nih.gov/15803924/)
91. van Gelder T et al. **Comparing mycophenolate mofetil regimens for de novo renal transplant recipients: the fixed-dose concentration-controlled trial.** *Transplantation* 2008. [PMID 18946341](https://pubmed.ncbi.nlm.nih.gov/18946341/)
92. Le Meur Y et al. **Individualized mycophenolate mofetil dosing based on drug exposure significantly improves patient outcomes after renal transplantation.** *Am J Transplant* 2007. [PMID 17908276](https://pubmed.ncbi.nlm.nih.gov/17908276/)
93. Sollinger HW **Mycophenolate mofetil for the prevention of acute rejection in primary cadaveric renal allograft recipients. U.S. Renal Transplant Mycophenolate Mofetil Study Group.** *Transplantation* 1995. [PMID 7645033](https://pubmed.ncbi.nlm.nih.gov/7645033/)
94. Sherwin CM et al. **The evolution of population pharmacokinetic models to describe the enterohepatic recycling of mycophenolic acid in solid organ transplantation and autoimmune disease.** *Clin Pharmacokinet* 2011. [PMID 21142265](https://pubmed.ncbi.nlm.nih.gov/21142265/)
95. Bergan S et al. **Personalized Therapy for Mycophenolate: Consensus Report by the International Association of Therapeutic Drug Monitoring and Clinical Toxicology.** *Ther Drug Monit* 2021. [PMID 33711005](https://pubmed.ncbi.nlm.nih.gov/33711005/)

## K. Corticosteroids · Early Withdrawal

> Corticosteroid PK/PD and early-withdrawal strategies.

96. Woodle ES et al. **Early Corticosteroid Cessation vs Long-term Corticosteroid Therapy in Kidney Transplant Recipients: Long-term Outcomes of a Randomized Clinical Trial.** *JAMA Surg* 2021. [PMID 33533901](https://pubmed.ncbi.nlm.nih.gov/33533901/)
97. Pascual J et al. **Steroid withdrawal in renal transplant patients on triple therapy with a calcineurin inhibitor and mycophenolate mofetil: a meta-analysis of randomized, controlled trials.** *Transplantation* 2004. [PMID 15599321](https://pubmed.ncbi.nlm.nih.gov/15599321/)
98. Haller MC et al. **Steroid avoidance or withdrawal for kidney transplant recipients.** *Cochrane Database Syst Rev* 2016. [PMID 27546100](https://pubmed.ncbi.nlm.nih.gov/27546100/)
99. Czock D et al. **Pharmacokinetics and pharmacodynamics of systemically administered glucocorticoids.** *Clin Pharmacokinet* 2005. [PMID 15634032](https://pubmed.ncbi.nlm.nih.gov/15634032/)

## L. Induction Therapy — ATG · basiliximab

> Induction therapy — polyclonal lymphocyte depletion (ATG) versus CD25 blockade (basiliximab).

100. Brennan DC et al. **Rabbit antithymocyte globulin versus basiliximab in renal transplantation.** *N Engl J Med* 2006. [PMID 17093248](https://pubmed.ncbi.nlm.nih.gov/17093248/)
101. Hill P et al. **Polyclonal and monoclonal antibodies for induction therapy in kidney transplant recipients.** *Cochrane Database Syst Rev* 2017. [PMID 28073178](https://pubmed.ncbi.nlm.nih.gov/28073178/)
102. Noël C et al. **Daclizumab versus antithymocyte globulin in high-immunological-risk renal transplant recipients.** *J Am Soc Nephrol* 2009. [PMID 19470677](https://pubmed.ncbi.nlm.nih.gov/19470677/)
103. Kho MM et al. **The effect of low and ultra-low dosages Thymoglobulin on peripheral T, B and NK cells in kidney transplant recipients.** *Transpl Immunol* 2012. [PMID 22410573](https://pubmed.ncbi.nlm.nih.gov/22410573/)
104. Höcker B et al. **Pharmacokinetics and immunodynamics of basiliximab in pediatric renal transplant recipients on mycophenolate mofetil comedication.** *Transplantation* 2008. [PMID 19005405](https://pubmed.ncbi.nlm.nih.gov/19005405/)

## M. Belatacept · Costimulation Blockade

> Belatacept and costimulation blockade, and resistance conferred by CD28-null memory T cells.

105. Larsen CP et al. **Rational development of LEA29Y (belatacept), a high-affinity variant of CTLA4-Ig with potent immunosuppressive properties.** *Am J Transplant* 2005. [PMID 15707398](https://pubmed.ncbi.nlm.nih.gov/15707398/)
106. Vincenti F et al. **A phase III study of belatacept-based immunosuppression regimens versus cyclosporine in renal transplant recipients (BENEFIT study).** *Am J Transplant* 2010. [PMID 20415897](https://pubmed.ncbi.nlm.nih.gov/20415897/)
107. Vincenti F et al. **Belatacept and Long-Term Outcomes in Kidney Transplantation.** *N Engl J Med* 2016. [PMID 26816011](https://pubmed.ncbi.nlm.nih.gov/26816011/)
108. Durrbach A et al. **Long-Term Outcomes in Belatacept- Versus Cyclosporine-Treated Recipients of Extended Criteria Donor Kidneys: Final Results From BENEFIT-EXT, a Phase III Randomized Study.** *Am J Transplant* 2016. [PMID 27130868](https://pubmed.ncbi.nlm.nih.gov/27130868/)
109. Bray RA et al. **De novo donor-specific antibodies in belatacept-treated vs cyclosporine-treated kidney-transplant recipients: Post hoc analyses of the randomized phase III BENEFIT and BENEFIT-EXT studies.** *Am J Transplant* 2018. [PMID 29509295](https://pubmed.ncbi.nlm.nih.gov/29509295/)
110. Espinosa J et al. **CD57(+) CD4 T Cells Underlie Belatacept-Resistant Allograft Rejection.** *Am J Transplant* 2016. [PMID 26603381](https://pubmed.ncbi.nlm.nih.gov/26603381/)

## N. mTOR Inhibitors · CNI Minimisation

> mTOR inhibitors and CNI-minimisation strategies.

111. Pascual J et al. **Everolimus with Reduced Calcineurin Inhibitor Exposure in Renal Transplantation.** *J Am Soc Nephrol* 2018. [PMID 29752413](https://pubmed.ncbi.nlm.nih.gov/29752413/)
112. Kahan BD **Efficacy of sirolimus compared with azathioprine for reduction of acute renal allograft rejection: a randomised multicentre study. The Rapamune US Study Group.** *Lancet* 2000. [PMID 10963197](https://pubmed.ncbi.nlm.nih.gov/10963197/)
113. Schena FP et al. **Conversion from calcineurin inhibitors to sirolimus maintenance therapy in renal allograft recipients: 24-month efficacy and safety results from the CONVERT trial.** *Transplantation* 2009. [PMID 19155978](https://pubmed.ncbi.nlm.nih.gov/19155978/)

## O. Treatment of Rejection · Desensitisation

> Treatment of rejection and desensitisation — PLEX/IVIG/rituximab/anti-IL-6/anti-CD38/imlifidase/bortezomib.

114. Sautenet B et al. **One-year Results of the Effects of Rituximab on Acute Antibody-Mediated Rejection in Renal Transplantation: RITUX ERAH, a Multicenter Double-blind Randomized Placebo-controlled Trial.** *Transplantation* 2016. [PMID 26555944](https://pubmed.ncbi.nlm.nih.gov/26555944/)
115. Choi J et al. **Assessment of Tocilizumab (Anti-Interleukin-6 Receptor Monoclonal) as a Potential Treatment for Chronic Antibody-Mediated Rejection and Transplant Glomerulopathy in HLA-Sensitized Renal Allograft Recipients.** *Am J Transplant* 2017. [PMID 28199785](https://pubmed.ncbi.nlm.nih.gov/28199785/)
116. Borski A et al. **Anti-interleukin-6 Antibody Clazakizumab in Antibody-mediated Renal Allograft Rejection: Accumulation of Antibody-neutralized Interleukin-6 Without Signs of Proinflammatory Rebound Phenomena.** *Transplantation* 2023. [PMID 35969004](https://pubmed.ncbi.nlm.nih.gov/35969004/)
117. Nickerson PW et al. **Clazakizumab for the treatment of chronic active antibody-mediated rejection (AMR) in kidney transplant recipients: Phase 3 IMAGINE study rationale and design.** *Trials* 2022. [PMID 36550562](https://pubmed.ncbi.nlm.nih.gov/36550562/)
118. Mayer KA et al. **Safety, tolerability, and efficacy of monoclonal CD38 antibody felzartamab in late antibody-mediated renal allograft rejection: study protocol for a phase 2 trial.** *Trials* 2022. [PMID 35395951](https://pubmed.ncbi.nlm.nih.gov/35395951/)
119. Jordan SC et al. **IgG Endopeptidase in Highly Sensitized Patients Undergoing Transplantation.** *N Engl J Med* 2017. [PMID 28767349](https://pubmed.ncbi.nlm.nih.gov/28767349/)
120. Kwun J et al. **Daratumumab in Sensitized Kidney Transplantation: Potentials and Limitations of Experimental and Clinical Use.** *J Am Soc Nephrol* 2019. [PMID 31227636](https://pubmed.ncbi.nlm.nih.gov/31227636/)
121. Montgomery RA et al. **Desensitization in HLA-incompatible kidney recipients and survival.** *N Engl J Med* 2011. [PMID 21793744](https://pubmed.ncbi.nlm.nih.gov/21793744/)
122. Orandi BJ et al. **Survival Benefit with Kidney Transplants from HLA-Incompatible Live Donors.** *N Engl J Med* 2016. [PMID 26962729](https://pubmed.ncbi.nlm.nih.gov/26962729/)
123. Velidedeoglu E et al. **Summary of 2017 FDA Public Workshop: Antibody-mediated Rejection in Kidney Transplantation.** *Transplantation* 2018. [PMID 29470345](https://pubmed.ncbi.nlm.nih.gov/29470345/)
124. Streichart L et al. **Tocilizumab in chronic active antibody-mediated rejection: rationale and protocol of an in-progress randomized controlled open-label multi-center trial (INTERCEPT study).** *Trials* 2024. [PMID 38519988](https://pubmed.ncbi.nlm.nih.gov/38519988/)

## P. BK Polyomavirus Nephropathy

> BK polyomavirus — the over-immunosuppression arm. The model's closed-loop BK screening controller.

125. Hirsch HH et al. **Prospective study of polyomavirus type BK replication and nephropathy in renal-transplant recipients.** *N Engl J Med* 2002. [PMID 12181403](https://pubmed.ncbi.nlm.nih.gov/12181403/)
126. Hirsch HH et al. **Polyomavirus-associated nephropathy in renal transplantation: interdisciplinary analyses and recommendations.** *Transplantation* 2005. [PMID 15912088](https://pubmed.ncbi.nlm.nih.gov/15912088/)
127. Hirsch HH et al. **BK polyomavirus in solid organ transplantation-Guidelines from the American Society of Transplantation Infectious Diseases Community of Practice.** *Clin Transplant* 2019. [PMID 30859620](https://pubmed.ncbi.nlm.nih.gov/30859620/)
128. Kant S et al. **BK Virus Nephropathy in Kidney Transplantation: A State-of-the-Art Review.** *Viruses* 2022. [PMID 35893681](https://pubmed.ncbi.nlm.nih.gov/35893681/)
129. Schaub S et al. **Reducing immunosuppression preserves allograft function in presumptive and definitive polyomavirus-associated nephropathy.** *Am J Transplant* 2010. [PMID 21114642](https://pubmed.ncbi.nlm.nih.gov/21114642/)
130. Kotton CN et al. **The Second International Consensus Guidelines on the Management of BK Polyomavirus in Kidney Transplantation.** *Transplantation* 2024. [PMID 38605438](https://pubmed.ncbi.nlm.nih.gov/38605438/)

## Q. CMV and Complications of Over-Immunosuppression

> CMV · other infections · malignancy · PTDM — the remaining price of over-immunosuppression.

131. Kotton CN et al. **The Third International Consensus Guidelines on the Management of Cytomegalovirus in Solid-organ Transplantation.** *Transplantation* 2018. [PMID 29596116](https://pubmed.ncbi.nlm.nih.gov/29596116/)
132. Humar A et al. **The efficacy and safety of 200 days valganciclovir cytomegalovirus prophylaxis in high-risk kidney transplant recipients.** *Am J Transplant* 2010. [PMID 20353469](https://pubmed.ncbi.nlm.nih.gov/20353469/)
133. Limaye AP et al. **Letermovir vs Valganciclovir for Prophylaxis of Cytomegalovirus in High-Risk Kidney Transplant Recipients: A Randomized Clinical Trial.** *JAMA* 2023. [PMID 37279999](https://pubmed.ncbi.nlm.nih.gov/37279999/)
134. Fishman JA **Infection in Organ Transplantation.** *Am J Transplant* 2017. [PMID 28117944](https://pubmed.ncbi.nlm.nih.gov/28117944/)
135. Sargen MR et al. **Spectrum of Nonkeratinocyte Skin Cancer Risk Among Solid Organ Transplant Recipients in the US.** *JAMA Dermatol* 2022. [PMID 35262623](https://pubmed.ncbi.nlm.nih.gov/35262623/)
136. Dharnidharka VR et al. **Post-transplant lymphoproliferative disorders.** *Nat Rev Dis Primers* 2016. [PMID 27189056](https://pubmed.ncbi.nlm.nih.gov/27189056/)
137. Sharif A et al. **Proceedings from an international consensus meeting on posttransplantation diabetes mellitus: recommendations and future directions.** *Am J Transplant* 2014. [PMID 25307034](https://pubmed.ncbi.nlm.nih.gov/25307034/)

## R. Medication Non-Adherence · Monitoring

> Medication non-adherence — the commonest correctable cause of late graft loss. The model's ADH state and its 'white-coat adherence' TDM assumption.

138. Schmid-Mohler G et al. **Non-adherence to immunosuppressive medication in renal transplant recipients within the scope of the Integrative Model of Behavioral Prediction: a cross-sectional study.** *Clin Transplant* 2010. [PMID 19674014](https://pubmed.ncbi.nlm.nih.gov/19674014/)
139. Butler JA et al. **Frequency and impact of nonadherence to immunosuppressants after renal transplantation: a systematic review.** *Transplantation* 2004. [PMID 15021846](https://pubmed.ncbi.nlm.nih.gov/15021846/)
140. Nevins TE et al. **Understanding Medication Nonadherence after Kidney Transplant.** *J Am Soc Nephrol* 2017. [PMID 28630231](https://pubmed.ncbi.nlm.nih.gov/28630231/)

## S. Diagnosis · Assessment of Renal Function · Guidelines

> Diagnosis · assessment of renal function · guidelines.

141. Inker LA et al. **New Creatinine- and Cystatin C-Based Equations to Estimate GFR without Race.** *N Engl J Med* 2021. [PMID 34554658](https://pubmed.ncbi.nlm.nih.gov/34554658/)
142. Kidney Disease: Improving Global Outcomes (KDIGO) Transplant Work Group **KDIGO clinical practice guideline for the care of kidney transplant recipients.** *Am J Transplant* 2009. [PMID 19845597](https://pubmed.ncbi.nlm.nih.gov/19845597/)
143. Chadban SJ et al. **KDIGO Clinical Practice Guideline on the Evaluation and Management of Candidates for Kidney Transplantation.** *Transplantation* 2020. [PMID 32301874](https://pubmed.ncbi.nlm.nih.gov/32301874/)
144. Rush D et al. **Beneficial effects of treatment of early subclinical rejection: a randomized study.** *J Am Soc Nephrol* 1998. [PMID 9808101](https://pubmed.ncbi.nlm.nih.gov/9808101/)
145. Hart A et al. **OPTN/SRTR 2015 Annual Data Report: Kidney.** *Am J Transplant* 2017. [PMID 28052609](https://pubmed.ncbi.nlm.nih.gov/28052609/)
146. Naesens M et al. **Surrogate Endpoints for Late Kidney Transplantation Failure.** *Transpl Int* 2022. [PMID 35669974](https://pubmed.ncbi.nlm.nih.gov/35669974/)

## T. Quantitative Systems Pharmacology · Modelling Methodology

> Quantitative systems pharmacology methodology and mrgsolve.

147. Rieger TR et al. **Improving the generation and selection of virtual populations in quantitative systems pharmacology models.** *Prog Biophys Mol Biol* 2018. [PMID 29902482](https://pubmed.ncbi.nlm.nih.gov/29902482/)
148. Cheng Y et al. **QSP Toolbox: Computational Implementation of Integrated Workflow Components for Deploying Multi-Scale Mechanistic Models.** *AAPS J* 2017. [PMID 28540623](https://pubmed.ncbi.nlm.nih.gov/28540623/)
149. Elmokadem A et al. **Quantitative Systems Pharmacology and Physiologically-Based Pharmacokinetic Modeling With mrgsolve: A Hands-On Tutorial.** *CPT Pharmacometrics Syst Pharmacol* 2019. [PMID 31652028](https://pubmed.ncbi.nlm.nih.gov/31652028/)
150. Nayak S et al. **Getting Innovative Therapies Faster to Patients at the Right Dose: Impact of Quantitative Pharmacology Towards First Registration and Expanding Therapeutic Use.** *Clin Pharmacol Ther* 2018. [PMID 29330855](https://pubmed.ncbi.nlm.nih.gov/29330855/)
151. Chan JR et al. **Role of Data in Development and Application of Quantitative Systems Pharmacology Models.** *Handb Exp Pharmacol* 2025. [PMID 40773021](https://pubmed.ncbi.nlm.nih.gov/40773021/)
152. Kirouac DC et al. **Reproducibility of Quantitative Systems Pharmacology Models: Current Challenges and Future Opportunities.** *CPT Pharmacometrics Syst Pharmacol* 2019. [PMID 30697975](https://pubmed.ncbi.nlm.nih.gov/30697975/)
153. Musante CJ et al. **Quantitative Systems Pharmacology: A Case for Disease Models.** *Clin Pharmacol Ther* 2017. [PMID 27709613](https://pubmed.ncbi.nlm.nih.gov/27709613/)

---

## Number of references

**153 papers** in total (PMID and title confirmed for every one of them by PubMed lookup).

## Quantitative targets used directly in model calibration

| Item | Target | Type of source |
|------|--------|-----------|
| 12-month eGFR (TAC/MPA/steroid) | 55-60 mL/min/1.73 m² | registration-trial control arms (Symphony, BENEFIT) |
| 12-month eGFR (belatacept, CNI-free) | 63-68 mL/min/1.73 m² | BENEFIT / BENEFIT-EXT |
| 12-month eGFR (cyclosporine) | 48-53 mL/min/1.73 m² | Symphony, BENEFIT control arms |
| eGFR slope in a stable graft | -1 ~ -2 mL/min/1.73 m²/year | cohorts |
| eGFR slope in chronic active ABMR | -5 ~ -10 mL/min/1.73 m²/year | cohorts |
| 5-year cumulative incidence of dnDSA | 15-25% (>40% with non-adherence) | Wiebe / Everly cohorts |
| Tacrolimus target trough | 8-10 at 0-3 months, 5-7 ng/mL thereafter | TDM consensus statements |
| Clearance in CYP3A5 expressers | roughly 1.5-2 times that of non-expressers | CPIC guideline |
| MPA AUC0-12 therapeutic window | 30-60 mg·h/L | TDM consensus statements |
| Tacrolimus IPV (trough CV%) | risk of dnDSA and graft loss rises above >30% | Borra / Shuker |
| BK viraemia → BKVAN threshold | sustained plasma 10⁴ copies/mL | AST / KDIGO guidance |
| dd-cfDNA cut-off | 0.5-1.0% (to rule out active rejection) | DART / Trifecta |
| Death-censored graft survival (deceased donor) | about 90% at 5 years, about 75% at 10 years | registries |

> ⚠️ This model is a **semi-quantitative QSP model for education and research**. The targets
> above are qualitative anchors for checking that the model lies in a clinically plausible
> range; they do not replace a formal fit to, and validation against, individual patient data.

