# Drug-Induced Liver Injury (DILI) — References
# Drug-Induced Liver Injury — Annotated Reference List

This list is the evidence for the structures, parameters, and clinical anchors used in `dili_qsp_model_en.dot` (the mechanistic map),
`dili_mrgsolve_model.R` (the ODE model), and `dili_shiny_app_en.R` (the dashboard).
The sections are divided so as to correspond to each subgraph cluster of the map and each ODE block of the model.

> **How to read this.** ★ marks a reference that determined the structure of the model itself (which terms exist);
> ▲ marks a reference that supplied a particular parameter value or a clinical calibration anchor.

---

## 1. Overview · epidemiology · causality assessment

1. ★ Andrade RJ, Chalasani N, Björnsson ES, et al. **Drug-induced liver injury.** *Nat Rev Dis Primers.* 2019;5(1):58. — the standard review of the whole pathophysiology of DILI. The 14 cluster compartments of this model follow the conceptual frame of this review. <https://pubmed.ncbi.nlm.nih.gov/31439850/>
2. ▲ Chalasani N, Bonkovsky HL, Fontana R, et al. **Features and Outcomes of 899 Patients With Drug-Induced Liver Injury: The DILIN Prospective Study.** *Gastroenterology.* 2015;148(7):1340-52.e7. — the anchor for the rate of chronicity (about 10–20%), the death/transplant rate, and the distribution of patterns. <https://pubmed.ncbi.nlm.nih.gov/25754159/>
3. ▲ Björnsson ES, Bergmann OM, Björnsson HK, et al. **Incidence, presentation, and outcomes in patients with drug-induced liver injury in the general population of Iceland.** *Gastroenterology.* 2013;144(7):1419-25. — a population-based incidence of 19.1/100,000/year. <https://pubmed.ncbi.nlm.nih.gov/23419359/>
4. ★ Danan G, Benichou C. **Causality assessment of adverse reactions to drugs—I. A novel method based on the conclusions of international consensus meetings: application to drug-induced liver injuries.** *J Clin Epidemiol.* 1993;46(11):1323-30. — the origin of RUCAM. <https://pubmed.ncbi.nlm.nih.gov/8229110/>
5. Hayashi PH, Lucena MI, Fontana RJ, et al. **A revised electronic version of RUCAM for the diagnosis of DILI.** *Hepatology.* 2022;76(1):18-31. — RECAM. <https://pubmed.ncbi.nlm.nih.gov/35014066/>
6. ★ Aithal GP, Watkins PB, Andrade RJ, et al. **Case definition and phenotype standardization in drug-induced liver injury.** *Clin Pharmacol Ther.* 2011;89(6):806-15. — the standard for the R ratio and for the definitions of the hepatocellular, cholestatic, and mixed patterns. The model's `R` output definition comes from here. <https://pubmed.ncbi.nlm.nih.gov/21544079/>
7. European Association for the Study of the Liver. **EASL Clinical Practice Guidelines: Drug-induced liver injury.** *J Hepatol.* 2019;70(6):1222-1261. <https://pubmed.ncbi.nlm.nih.gov/30926241/>
8. Fontana RJ, Liou I, Reuben A, et al. **AASLD practice guidance on drug, herbal, and dietary supplement-induced liver injury.** *Hepatology.* 2023;77(3):1036-1065. <https://pubmed.ncbi.nlm.nih.gov/35899384/>
9. ▲ Navarro VJ, Khan I, Björnsson E, et al. **Liver injury from herbal and dietary supplements.** *Hepatology.* 2017;65(1):363-373. — HDS accounts for more than 20% of the US DILIN registry. <https://pubmed.ncbi.nlm.nih.gov/27677775/>

---

## 2. Hy's Law · prognosis · the regulatory framework

10. ★▲ Temple R. **Hy's law: predicting serious hepatotoxicity.** *Pharmacoepidemiol Drug Saf.* 2006;15(4):241-3. — ALT ≥3×ULN **AND** TBIL ≥2×ULN → death/transplant about 10%. The rule that the central claim of this model ("Hy's Law is a joint test of rate × reserve") is aimed at. <https://pubmed.ncbi.nlm.nih.gov/16552790/>
11. ★ Robles-Diaz M, Lucena MI, Kaplowitz N, et al. **Use of Hy's law and a new composite algorithm to predict acute liver failure in patients with drug-induced liver injury.** *Gastroenterology.* 2014;147(1):109-118.e5. — the predictive performance of Hy's Law and its improvement on an nR basis. <https://pubmed.ncbi.nlm.nih.gov/24704526/>
12. ▲ U.S. FDA. **Guidance for Industry: Drug-Induced Liver Injury — Premarketing Clinical Evaluation.** 2009. — the regulatory case definition of Hy's Law and the stopping rule. <https://www.fda.gov/regulatory-information/search-fda-guidance-documents/drug-induced-liver-injury-premarketing-clinical-evaluation>
13. ▲ O'Grady JG, Alexander GJ, Hayllar KM, Williams R. **Early indicators of prognosis in fulminant hepatic failure.** *Gastroenterology.* 1989;97(2):439-45. — the origin of the King's College criteria (APAP form: pH<7.30, or INR>6.5 + Cr>3.4 mg/dL + grade III–IV encephalopathy). <https://pubmed.ncbi.nlm.nih.gov/2490426/>
14. Bernal W, Wendon J. **Acute liver failure.** *N Engl J Med.* 2013;369(26):2525-34. <https://pubmed.ncbi.nlm.nih.gov/24369077/>
15. ▲ Lee WM. **Acetaminophen (APAP) hepatotoxicity—Isn't it time for APAP to go away?** *J Hepatol.* 2017;67(6):1324-1331. — the largest cause of ALF in the United States (about 46%). <https://pubmed.ncbi.nlm.nih.gov/28734939/>

---

## 3. Acetaminophen biotransformation and dose thresholds

16. ★▲ Mitchell JR, Jollow DJ, Potter WZ, Gillette JR, Brodie BB. **Acetaminophen-induced hepatic necrosis. IV. Protective role of glutathione.** *J Pharmacol Exp Ther.* 1973;187(1):211-7. — the classic paper establishing that GSH depletion is a precondition of necrosis. The basis for the "GSH gate" structure of the model (`v_bind` competitively inhibited by GSH). <https://pubmed.ncbi.nlm.nih.gov/4746329/>
17. ★ Jollow DJ, Mitchell JR, Potter WZ, et al. **Acetaminophen-induced hepatic necrosis. II. Role of covalent binding in vivo.** *J Pharmacol Exp Ther.* 1973;187(1):195-202. <https://pubmed.ncbi.nlm.nih.gov/4746327/>
18. ★▲ Dahlin DC, Miwa GT, Lu AY, Nelson SD. **N-acetyl-p-benzoquinone imine: a cytochrome P-450-mediated oxidation product of acetaminophen.** *Proc Natl Acad Sci USA.* 1984;81(5):1327-31. — identification of NAPQI. <https://pubmed.ncbi.nlm.nih.gov/6424115/>
19. ▲ Prescott LF. **Kinetics and metabolism of paracetamol and phenacetin.** *Br J Clin Pharmacol.* 1980;10(Suppl 2):291S-298S. — at therapeutic doses, glucuronide ≈ 50–60%, sulfate ≈ 25–35%, the CYP route ≈ 5–10%, renal excretion ≈ 5%. The anchor for the Vmax/Km partitioning of the model and for its total clearance (about 19–21 L/h). <https://pubmed.ncbi.nlm.nih.gov/7002186/>
20. ★▲ Slattery JT, Wilson JM, Kalhorn TF, Nelson SD. **Dose-dependent pharmacokinetics of acetaminophen: evidence of glutathione depletion in humans.** *Clin Pharmacol Ther.* 1987;41(4):413-8. — dose-dependent saturation of sulfation in humans. The **PAPS/UDPGA cofactor depletion ODEs** of the model come from here (this term is what makes "at a high dose the activated fraction rises of its own accord" come out of the calculation). <https://pubmed.ncbi.nlm.nih.gov/3829578/>
21. ▲ Rumack BH, Matthew H. **Acetaminophen poisoning and toxicity.** *Pediatrics.* 1975;55(6):871-6. — the origin of the Rumack-Matthew nomogram (the 4-hour 200 µg/mL line; in the United States a 150 µg/mL 'treatment line'). <https://pubmed.ncbi.nlm.nih.gov/1134886/>
22. ▲ Rumack BH. **Acetaminophen hepatotoxicity: the first 35 years.** *J Toxicol Clin Toxicol.* 2002;40(1):3-20. — the 150 mg/kg / 7.5 g threshold, with risk rising steeply above 250 mg/kg. The clinical anchor for the position of the knee in the model's dose-response. <https://pubmed.ncbi.nlm.nih.gov/11990202/>
23. ▲ Zhao L, Pickering G. **Paracetamol metabolism and related genetic differences.** *Drug Metab Rev.* 2011;43(1):41-52. <https://pubmed.ncbi.nlm.nih.gov/20977384/>
24. Chiew AL, Reith D, Pomerleau A, et al. **Updated guidelines for the management of paracetamol poisoning in Australia and New Zealand.** *Med J Aust.* 2020;212(4):175-183. <https://pubmed.ncbi.nlm.nih.gov/31786822/>

---

## 4. Glutathione, cysteine, thiol defence

25. ★ Lu SC. **Glutathione synthesis.** *Biochim Biophys Acta.* 2013;1830(5):3143-53. — GCL is the rate-limiting enzyme and **cysteine the rate-limiting substrate**. Why the effect of NAC enters the model only through the `CYS` state variable. <https://pubmed.ncbi.nlm.nih.gov/22995213/>
26. ★▲ Lu SC. **Regulation of hepatic glutathione synthesis: current concepts and controversies.** *FASEB J.* 1999;13(10):1169-83. — a hepatic GSH concentration of 5–10 mM and a turnover half-life of 2–4 hours. <https://pubmed.ncbi.nlm.nih.gov/10385608/>
27. ★ Fernández-Checa JC, Kaplowitz N. **Hepatic mitochondrial glutathione: transport and role in disease and toxicity.** *Toxicol Appl Pharmacol.* 2005;204(3):263-73. — the mitochondrial GSH pool is a separate compartment and its depletion is lethal. <https://pubmed.ncbi.nlm.nih.gov/15845418/>
28. ▲ Chen Y, Dong H, Thompson DC, et al. **Glutathione defense mechanism in liver injury: insights from animal models.** *Food Chem Toxicol.* 2013;60:38-44. <https://pubmed.ncbi.nlm.nih.gov/23856494/>
29. ★ Yuan L, Kaplowitz N. **Glutathione in liver diseases and hepatotoxicity.** *Mol Aspects Med.* 2009;30(1-2):29-41. <https://pubmed.ncbi.nlm.nih.gov/18786561/>

---

## 5. Mitochondrial toxicity and bioenergetics

30. ★ Jaeschke H, McGill MR, Ramachandran A. **Oxidant stress, mitochondria, and cell death mechanisms in drug-induced liver injury: lessons learned from acetaminophen hepatotoxicity.** *Drug Metab Rev.* 2012;44(1):88-106. — the standard account of the mitochondrial adduct → ROS → MPT → necrosis pathway. The skeleton of cluster 4 of the model. <https://pubmed.ncbi.nlm.nih.gov/22229890/>
31. ★▲ Kon K, Kim JS, Jaeschke H, Lemasters JJ. **Mitochondrial permeability transition in acetaminophen-induced necrosis and apoptosis of cultured mouse hepatocytes.** *Hepatology.* 2004;40(5):1170-9. — **with ATP, apoptosis; without it, necrosis**. The basis for the ATP gate `1/(1+(ATP/K)^4)` in the MPT term of the model. <https://pubmed.ncbi.nlm.nih.gov/15486922/>
32. ★ Ramachandran A, Jaeschke H. **Acetaminophen Hepatotoxicity.** *Semin Liver Dis.* 2019;39(2):221-234. <https://pubmed.ncbi.nlm.nih.gov/30849782/>
33. ★ Pessayre D, Fromenty B, Berson A, et al. **Central role of mitochondria in drug-induced liver injury.** *Drug Metab Rev.* 2012;44(1):34-87. — the three independent axes of mitochondrial injury: uncoupling, inhibition of β-oxidation, and mtDNA depletion. The model's `FUNC` (uncoupling/FAO-blockade burden) parameter. <https://pubmed.ncbi.nlm.nih.gov/21892896/>
34. ▲ Nadanaciva S, Will Y. **New insights in drug-induced mitochondrial toxicity.** *Curr Pharm Des.* 2011;17(20):2100-12. <https://pubmed.ncbi.nlm.nih.gov/21718246/>
35. ★ Masubuchi Y, Suda C, Horie T. **Involvement of mitochondrial permeability transition in acetaminophen-induced liver injury in mice.** *J Hepatol.* 2005;42(1):110-6. <https://pubmed.ncbi.nlm.nih.gov/15629515/>

---

## 6. The JNK–Sab amplification loop — the source of bistability

36. ★★ Hanawa N, Shinohara M, Saberi B, et al. **Role of JNK translocation to mitochondria leading to inhibition of mitochondria bioenergetics in acetaminophen-induced liver injury.** *J Biol Chem.* 2008;283(20):13565-77. — **p-JNK translocates to the mitochondrion, inhibits respiration, and makes still more ROS** = the positive feedback loop of this model itself. <https://pubmed.ncbi.nlm.nih.gov/18337250/>
37. ★★ Win S, Than TA, Han D, Petrovic LM, Kaplowitz N. **c-Jun N-terminal kinase (JNK)-dependent acute liver injury from acetaminophen or tumor necrosis factor (TNF) requires mitochondrial Sab protein expression in mice.** *J Biol Chem.* 2011;286(40):35071-8. — proves that Sab (SH3BP5) is the docking protein of the loop. Delete Sab and the injury disappears → the **necessity** of the loop. <https://pubmed.ncbi.nlm.nih.gov/21844199/>
38. ★★ Win S, Than TA, Zhang J, Oo C, Min RWM, Kaplowitz N. **New insights into the role and mechanism of c-Jun-N-terminal kinase signaling in the pathobiology of liver diseases.** *Hepatology.* 2018;67(5):2013-2024. — transient JNK activation is a survival signal, **sustained** activation a death signal. The clinical reality the model expresses as "bistable". <https://pubmed.ncbi.nlm.nih.gov/29194686/>
39. ★ Gunawan BK, Liu ZX, Han D, Hanawa N, Gaarde WA, Kaplowitz N. **c-Jun N-terminal kinase plays a major role in murine acetaminophen hepatotoxicity.** *Gastroenterology.* 2006;131(1):165-78. <https://pubmed.ncbi.nlm.nih.gov/16831600/>
40. ★ Nakagawa H, Maeda S, Hikiba Y, et al. **Deletion of apoptosis signal-regulating kinase 1 attenuates acetaminophen-induced liver injury by inhibiting c-Jun N-terminal kinase activation.** *Gastroenterology.* 2008;135(4):1311-21. — ASK1 is the ROS→JNK relay. <https://pubmed.ncbi.nlm.nih.gov/18700144/>
41. Saito C, Lemasters JJ, Jaeschke H. **c-Jun N-terminal kinase modulates oxidant stress and peroxynitrite formation independent of inducible nitric oxide synthase in acetaminophen hepatotoxicity.** *Toxicol Appl Pharmacol.* 2010;246(1-2):8-17. <https://pubmed.ncbi.nlm.nih.gov/20423716/>
42. ★ Win S, Min RWM, Zhang J, et al. **Hepatic Mitochondrial SAB Deletion or Knockdown Alleviates Diet-Induced Metabolic Syndrome, Steatohepatitis, and Hepatic Fibrosis.** *Hepatology.* 2021;74(6):3127-3145. — the loop is common to liver diseases beyond APAP. <https://pubmed.ncbi.nlm.nih.gov/34272738/>

---

## 7. Nrf2 adaptation and autophagy

43. ★ Enomoto A, Itoh K, Nagayoshi E, et al. **High sensitivity of Nrf2 knockout mice to acetaminophen hepatotoxicity associated with decreased expression of ARE-regulated drug metabolizing enzymes and antioxidant genes.** *Toxicol Sci.* 2001;59(1):169-77. — susceptibility rises steeply in Nrf2-null mice → the basis for `NRF2` in the model raising both GSH synthetic capacity and ROS scavenging capacity multiplicatively. <https://pubmed.ncbi.nlm.nih.gov/11134556/>
44. ★ Goldring CE, Kitteringham NR, Elsby R, et al. **Activation of hepatic Nrf2 in vivo by acetaminophen in CD-1 mice.** *Hepatology.* 2004;39(5):1267-76. <https://pubmed.ncbi.nlm.nih.gov/15122754/>
45. ★ Ni HM, Bockus A, Boggess N, Jaeschke H, Ding WX. **Activation of autophagy protects against acetaminophen-induced hepatotoxicity.** *Hepatology.* 2012;55(1):222-32. — autophagy clears damaged mitochondria and adducts. The model's `KADD_REP × autophagy(ATP)` term. <https://pubmed.ncbi.nlm.nih.gov/21932416/>
46. ▲ Dinkova-Kostova AT, Kostov RV, Kazantsev AG. **The role of Nrf2 signaling in counteracting neurodegenerative diseases.** *FEBS J.* 2018;285(19):3576-3590. — a general account of the KEAP1 cysteine sensor mechanism. <https://pubmed.ncbi.nlm.nih.gov/29323772/>
47. ★ Kaplowitz N. **Adaptation vs. injury: what determines the outcome?** In: *Drug-Induced Liver Disease* (3rd ed). — the concept that adaptation is the outcome for the great majority of those exposed. (Overview: Watkins PB. **Idiosyncratic liver injury: challenges and approaches.** *Toxicol Pathol.* 2005;33(1):1-5. <https://pubmed.ncbi.nlm.nih.gov/15805049/>)

---

## 8. Bile acids, BSEP, cholestatic DILI

48. ★★ Morgan RE, Trauner M, van Staden CJ, et al. **Interference with bile salt export pump function is a susceptibility factor for human liver injury in drug development.** *Toxicol Sci.* 2010;118(2):485-500. — the lower the BSEP inhibition IC50, the higher the DILI risk. The model's `KI_BSEP` competitive inhibition term. <https://pubmed.ncbi.nlm.nih.gov/20829430/>
49. ★▲ Dawson S, Stahl S, Paul N, Barber J, Kenna JG. **In vitro inhibition of the bile salt export pump correlates with risk of cholestatic drug-induced liver injury in humans.** *Drug Metab Dispos.* 2012;40(1):130-8. — a raised risk of cholestatic DILI for drugs with a BSEP IC50 < 25–50 µM. <https://pubmed.ncbi.nlm.nih.gov/21965623/>
50. ★ Stieger B, Fattinger K, Madon J, Kullak-Ublick GA, Meier PJ. **Drug- and estrogen-induced cholestasis through inhibition of the hepatocellular bile salt export pump (Bsep) of rat liver.** *Gastroenterology.* 2000;118(2):422-30. <https://pubmed.ncbi.nlm.nih.gov/10648470/>
51. ★ Woolbright BL, Jaeschke H. **Novel insight into mechanisms of cholestatic liver injury.** *World J Gastroenterol.* 2012;18(36):4985-93. — the cytotoxicity of bile acids is not simple detergent action but is inflammation-mediated (induction of CXCL1/2). The `BADETER → CXCL` edge in the model. <https://pubmed.ncbi.nlm.nih.gov/23049205/>
52. ★ Chiang JYL. **Bile acid metabolism and signaling.** *Compr Physiol.* 2013;3(3):1191-212. — the FXR-SHP-CYP7A1 negative feedback and the intestinal FGF19 axis. The model's `FMAX_FXR` partial inhibition. <https://pubmed.ncbi.nlm.nih.gov/23897684/>
53. ▲ Slijepcevic D, van de Graaf SFJ. **Bile Acid Uptake Transporters as Targets for Therapy.** *Dig Dis.* 2017;35(3):251-258. <https://pubmed.ncbi.nlm.nih.gov/28249287/>
54. ★ Aleo MD, Luo Y, Swiss R, Bonin PD, Potter DM, Will Y. **Human drug-induced liver injury severity is highly associated with dual inhibition of liver mitochondrial function and bile salt export pump.** *Hepatology.* 2014;60(3):1015-22. — **the double hit of mitochondrial toxicity + BSEP inhibition** is the worst combination. Why BSEP activity is scaled `× ATP` in the model (hepatocellular injury creates a secondary cholestasis). <https://pubmed.ncbi.nlm.nih.gov/24799086/>
55. Padda MS, Sanchez M, Akhtar AJ, Boyer JL. **Drug-induced cholestasis.** *Hepatology.* 2011;53(4):1377-87. <https://pubmed.ncbi.nlm.nih.gov/21480339/>

---

## 9. Innate immunity and sterile inflammation

56. ★ Jaeschke H, Ramachandran A. **Mechanisms and pathophysiological significance of sterile inflammation during acetaminophen hepatotoxicity.** *Food Chem Toxicol.* 2020;138:111240. — the DAMP → Kupffer cell → cytokine axis both **amplifies the injury and initiates the repair** (the double edge of TNF). <https://pubmed.ncbi.nlm.nih.gov/32145352/>
57. ★ Scaffidi P, Misteli T, Bianchi ME. **Release of chromatin protein HMGB1 by necrotic cells triggers inflammation.** *Nature.* 2002;418(6894):191-5. — HMGB1 = the representative DAMP. <https://pubmed.ncbi.nlm.nih.gov/12110890/>
58. ★ Imaeda AB, Watanabe A, Sohail MA, et al. **Acetaminophen-induced hepatotoxicity in mice is dependent on Tlr9 and the Nalp3 inflammasome.** *J Clin Invest.* 2009;119(2):305-14. <https://pubmed.ncbi.nlm.nih.gov/19164858/>
59. ★ Bourdi M, Masubuchi Y, Reilly TP, et al. **Protection against acetaminophen-induced liver injury and lethality by interleukin 10: role of inducible nitric oxide synthase.** *Hepatology.* 2002;35(2):289-98. — the protective effect of IL-10. The model's `IL10 ⊣ KC, TNF` inhibitory edges. <https://pubmed.ncbi.nlm.nih.gov/11826401/>
60. ★ Yamada Y, Kirillova I, Peschon JJ, Fausto N. **Initiation of liver growth by tumor necrosis factor: deficient liver regeneration in mice lacking type I tumor necrosis factor receptor.** *Proc Natl Acad Sci USA.* 1997;94(4):1441-6. — TNF/IL-6 prime regeneration. The model's `prim` term. <https://pubmed.ncbi.nlm.nih.gov/9037072/>
61. Antoniades CG, Quaglia A, Taams LS, et al. **Source and characterization of hepatic macrophages in acetaminophen-induced acute liver failure in humans.** *Hepatology.* 2012;56(2):735-46. <https://pubmed.ncbi.nlm.nih.gov/22334567/>
62. Woolbright BL, Jaeschke H. **Role of the inflammasome in acetaminophen-induced liver injury and acute liver failure.** *J Hepatol.* 2017;66(4):836-848. <https://pubmed.ncbi.nlm.nih.gov/27913221/>

---

## 10. Adaptive immunity, HLA, idiosyncratic DILI

63. ★★ Uetrecht J, Naisbitt DJ. **Idiosyncratic adverse drug reactions: current concepts.** *Pharmacol Rev.* 2013;65(2):779-808. — the hapten hypothesis, the p-i concept, and the requirement for a danger signal. Why the T-cell activation term of the model is the product of **adduct × danger signal × (1/tolerance)**. <https://pubmed.ncbi.nlm.nih.gov/23476052/>
64. ★★ Daly AK, Donaldson PT, Bhatnagar P, et al. **HLA-B*5701 genotype is a major determinant of drug-induced liver injury due to flucloxacillin.** *Nat Genet.* 2009;41(7):816-9. — an odds ratio of about 80. The model's `HLA` switch. <https://pubmed.ncbi.nlm.nih.gov/19483685/>
65. ★▲ Lucena MI, Molokhia M, Shen Y, et al. **Susceptibility to amoxicillin-clavulanate-induced liver injury is influenced by multiple HLA class I and II alleles.** *Gastroenterology.* 2011;141(1):338-47. — HLA-DRB1*15:01 and others. <https://pubmed.ncbi.nlm.nih.gov/21570397/>
66. ★ Urban TJ, Nicoletti P, Chalasani N, et al. **Minocycline hepatotoxicity: Clinical characterization and identification of HLA-B*35:02 as a risk factor.** *J Hepatol.* 2017;67(1):137-144. <https://pubmed.ncbi.nlm.nih.gov/28323125/>
67. ★★ Metushi IG, Hayes MA, Uetrecht J. **Treatment of PD-1−/− mice with anti-CTLA-4 antibody unmasks liver injury from amodiaquine.** *Hepatology.* 2015;61(4):1332-42. — direct proof that **removing immune tolerance converts a latent adduct burden into clinical hepatitis**. The structure of the ICI scenario of the model comes from here. <https://pubmed.ncbi.nlm.nih.gov/25482010/>
68. ★ Chakraborty M, Fullerton AM, Semple K, et al. **Drug-induced allergic hepatitis develops in mice when myeloid-derived suppressor cells are depleted prior to halothane treatment.** *Hepatology.* 2015;62(2):546-57. <https://pubmed.ncbi.nlm.nih.gov/25712247/>
68b. Cho T, Uetrecht J. **How Reactive Metabolites Induce an Immune Response That Sometimes Leads to an Idiosyncratic Drug Reaction.** *Chem Res Toxicol.* 2017;30(1):295-314. <https://pubmed.ncbi.nlm.nih.gov/27775332/>
69. ▲ De Martin E, Michot JM, Papouin B, et al. **Characterization of liver injury induced by cancer immunotherapy using immune checkpoint inhibitors.** *J Hepatol.* 2018;68(6):1181-1190. — the clinical and pathological phenotype of ICI hepatitis and its steroid responsiveness. <https://pubmed.ncbi.nlm.nih.gov/29427729/>
70. ▲ Peeraphatdit TB, Wang J, Odenwald MA, Hu S, Hart J, Charlton MR. **Hepatotoxicity From Immune Checkpoint Inhibitors: A Systematic Review and Management Recommendation.** *Hepatology.* 2020;72(1):315-329. — the incidence of grade 3–4 hepatotoxicity and stepped steroid/MMF therapy. <https://pubmed.ncbi.nlm.nih.gov/32167613/>

---

## 11. Mechanistic biomarkers

71. ★★ Antoine DJ, Dear JW, Lewis PS, et al. **Mechanistic biomarkers provide early and sensitive detection of acetaminophen-induced acute liver injury at first presentation to hospital.** *Hepatology.* 2013;58(2):777-87. — miR-122, HMGB1, and K18 rise before ALT. Because of the short half-life of miR-122, **that earlier rise comes out of the calculation** in the model. <https://pubmed.ncbi.nlm.nih.gov/23390034/>
72. ★▲ Starkey Lewis PJ, Dear J, Platt V, et al. **Circulating microRNAs as potential markers of human drug-induced liver injury.** *Hepatology.* 2011;54(5):1767-76. <https://pubmed.ncbi.nlm.nih.gov/22045675/>
73. ▲ Church RJ, Kullak-Ublick GA, Aubrecht J, et al. **Candidate biomarkers for the diagnosis and prognosis of drug-induced liver injury: An international collaborative effort.** *Hepatology.* 2019;69(2):760-773. — multicentre validation. <https://pubmed.ncbi.nlm.nih.gov/29357190/>
74. ★▲ McGill MR, Sharpe MR, Williams CD, Taha M, Curry SC, Jaeschke H. **The mechanism underlying acetaminophen-induced hepatotoxicity in humans and mice involves mitochondrial damage and nuclear DNA fragmentation.** *J Clin Invest.* 2012;122(4):1574-83. — the rise in GLDH, mtDNA, and nDNA fragments in humans. <https://pubmed.ncbi.nlm.nih.gov/22378043/>
75. ▲ Dear JW, Clarke JI, Francis B, et al. **Risk stratification after paracetamol overdose using mechanistic biomarkers: results from two prospective cohort studies.** *Lancet Gastroenterol Hepatol.* 2018;3(2):104-113. <https://pubmed.ncbi.nlm.nih.gov/29146439/>
76. ▲ James LP, Letzig L, Simpson PM, et al. **Pharmacokinetics of acetaminophen-protein adducts in adults with acetaminophen overdose and acute liver failure.** *Drug Metab Dispos.* 2009;37(8):1779-84. — the kinetics of circulating APAP-Cys adducts (a confirmatory biomarker of exposure). <https://pubmed.ncbi.nlm.nih.gov/19439490/>
77. ▲ Ozer J, Ratner M, Shaw M, Bailey W, Schomaker S. **The current state of serum biomarkers of hepatotoxicity.** *Toxicology.* 2008;245(3):194-205. — a serum half-life of about 47 hours for ALT and about 17 hours for AST. The model's `KALT_EL` and `KAST_EL`. <https://pubmed.ncbi.nlm.nih.gov/18291570/>

---

## 12. N-acetylcysteine and other therapeutics

78. ★★ Prescott LF, Illingworth RN, Critchley JA, Stewart MJ, Adam RD, Proudfoot AT. **Intravenous N-acetylcystine: the treatment of choice for paracetamol poisoning.** *Br Med J.* 1979;2(6198):1097-100. — the origin of the Prescott regimen (150 → 50 → 100 mg/kg). The NAC infusion schedule of the model. <https://pubmed.ncbi.nlm.nih.gov/519312/>
79. ★★ Smilkstein MJ, Knapp GL, Kulig KW, Rumack BH. **Efficacy of oral N-acetylcysteine in the treatment of acetaminophen overdose. Analysis of the national multicenter study (1976 to 1985).** *N Engl J Med.* 1988;319(24):1557-62. — **given within 10 hours, hepatotoxicity is close to 0%; between 16 and 24 hours it deteriorates sharply.** The central clinical anchor that the NAC time-window simulation of this model has to reproduce. <https://pubmed.ncbi.nlm.nih.gov/3059186/>
80. ▲ Bateman DN, Dear JW, Thanacoody HK, et al. **Reduction of adverse effects from intravenous acetylcysteine treatment for paracetamol poisoning: a randomised controlled trial.** *Lancet.* 2014;383(9918):697-704. — the SNAP 12-hour regimen. <https://pubmed.ncbi.nlm.nih.gov/24290406/>
81. ★ Lauterburg BH, Corcoran GB, Mitchell JR. **Mechanism of action of N-acetylcysteine in the protection against the hepatotoxicity of acetaminophen in rats in vivo.** *J Clin Invest.* 1983;71(4):980-91. — the principal mechanism of NAC is **resynthesis of GSH through the supply of cysteine**. <https://pubmed.ncbi.nlm.nih.gov/6833497/>
82. ★ Lee WM, Hynan LS, Rossaro L, et al. **Intravenous N-acetylcysteine improves transplant-free survival in early stage non-acetaminophen acute liver failure.** *Gastroenterology.* 2009;137(3):856-64. — beneficial in non-APAP ALF as well (a microcirculatory/oxygen delivery mechanism). <https://pubmed.ncbi.nlm.nih.gov/19524577/>
83. ▲ Akakpo JY, Ramachandran A, Kandel SE, et al. **4-Methylpyrazole protects against acetaminophen hepatotoxicity in mice and in primary human hepatocytes.** *Hum Exp Toxicol.* 2018;37(12):1310-1322. — fomepizole (CYP2E1 + JNK inhibition). <https://pubmed.ncbi.nlm.nih.gov/29768939/>
84. ▲ Morrison EE, Oatey K, Gallagher B, et al. **Principal results of a randomised open label exploratory, safety and tolerability study with calmangafodipir in patients treated with a 12 h regimen of N-acetylcysteine for paracetamol overdose (POP trial).** *EBioMedicine.* 2019;46:423-430. <https://pubmed.ncbi.nlm.nih.gov/31351929/>
85. ▲ Larsen FS, Schmidt LE, Bernsmeier C, et al. **High-volume plasma exchange in patients with acute liver failure: An open randomised controlled trial.** *J Hepatol.* 2016;64(1):69-78. — high-volume plasma exchange improves survival in ALF. <https://pubmed.ncbi.nlm.nih.gov/26325537/>
86. ▲ Lheureux PE, Penaloza A, Zahir S, Gris M. **Science review: carnitine in the treatment of valproic acid-induced toxicity.** *Crit Care.* 2005;9(5):431-40. <https://pubmed.ncbi.nlm.nih.gov/16277730/>

---

## 13. Host risk factors and physicochemical rules

87. ★▲ Chen M, Borlak J, Tong W. **High lipophilicity and high daily dose of oral medications are associated with significant risk for drug-induced liver injury.** *Hepatology.* 2013;58(1):388-96. — the "Rule-of-Two" (logP > 3 and a daily dose > 100 mg). The `LIPOPHIL` node of the model's map. <https://pubmed.ncbi.nlm.nih.gov/23258593/>
88. ★▲ Lammert C, Einarsson S, Saha C, Niklasson A, Bjornsson E, Chalasani N. **Relationship between daily dose of oral medications and idiosyncratic drug-induced liver injury: search for signals.** *Hepatology.* 2008;47(6):2003-9. — the risk concentrates at daily doses ≥50 mg. <https://pubmed.ncbi.nlm.nih.gov/18454504/>
89. ★ Whitcomb DC, Block GD. **Association of acetaminophen hepatotoxicity with fasting and ethanol use.** *JAMA.* 1994;272(23):1845-50. — the double hit of fasting and alcohol. The vulnerable-host scenario of the model (`FCYP`↑, `FGSH`↓, `CYSBASE`↓). <https://pubmed.ncbi.nlm.nih.gov/7990219/>
90. ★ Zimmerman HJ, Maddrey WC. **Acetaminophen (paracetamol) hepatotoxicity with regular intake of alcohol: analysis of instances of therapeutic misadventure.** *Hepatology.* 1995;22(3):767-73. <https://pubmed.ncbi.nlm.nih.gov/7657281/>
91. ▲ Michaut A, Moreau C, Robin MA, Fromenty B. **Acetaminophen-induced liver injury in obesity and nonalcoholic fatty liver disease.** *Liver Int.* 2014;34(7):e171-9. — the raised CYP2E1 and reduced mitochondrial reserve of obesity/MASLD. <https://pubmed.ncbi.nlm.nih.gov/24575897/>
92. ▲ Nicoletti P, Aithal GP, Bjornsson ES, et al. **Association of Liver Injury From Specific Drugs, or Groups of Drugs, With Polymorphisms in HLA and Other Genes in a Genome-Wide Association Study.** *Gastroenterology.* 2017;152(5):1078-1089. — non-HLA factors such as PTPN22. <https://pubmed.ncbi.nlm.nih.gov/28043905/>

---

## 14. Regeneration, functional reserve, chronicity

93. ★ Michalopoulos GK, Bhushan B. **Liver regeneration: biological and pathological mechanisms and implications.** *Nat Rev Gastroenterol Hepatol.* 2021;18(1):40-55. — HGF/c-Met-driven regeneration and the termination signal of TGF-β1. The structure of the model's `regen` term. <https://pubmed.ncbi.nlm.nih.gov/32764740/>
94. ★ Fausto N, Campbell JS, Riehle KJ. **Liver regeneration.** *Hepatology.* 2006;43(2 Suppl 1):S45-53. <https://pubmed.ncbi.nlm.nih.gov/16447274/>
95. ▲ Fontana RJ, Hayashi PH, Barnhart H, et al. **Persistent liver biochemistry abnormalities are more common in older patients and those with cholestatic drug induced liver injury.** *Am J Gastroenterol.* 2015;110(10):1450-9. — chronicity is common in the cholestatic pattern and in older patients. <https://pubmed.ncbi.nlm.nih.gov/26346867/>
96. ▲ Bonkovsky HL, Kleiner DE, Gu J, et al. **Clinical presentations and outcomes of bile duct loss caused by drugs and herbal and dietary supplements.** *Hepatology.* 2017;65(4):1267-1277. — vanishing bile duct syndrome (VBDS). <https://pubmed.ncbi.nlm.nih.gov/27981596/>

---

## 15. QSP modelling methodology · DILIsym and prior models

97. ★★ Howell BA, Yang Y, Kumar R, et al. **In vitro to in vivo extrapolation and species response comparisons for drug-induced liver injury (DILI) using DILIsym™: a mechanistic, mathematical model of DILI.** *J Pharmacokinet Pharmacodyn.* 2012;39(5):527-41. — the standard platform for DILI QSP. This model refers to the published structural concepts of DILIsym (the GSH, mitochondrial, bile acid, and life-cycle submodels) but simplifies and reimplements them independently. <https://pubmed.ncbi.nlm.nih.gov/22996471/>
98. ★★ Watkins PB. **DILIsym: Quantitative systems toxicology impacting drug development.** *Curr Opin Toxicol.* 2020;23-24:67-73. <https://doi.org/10.1016/j.cotox.2020.06.003>
99. ★ Woodhead JL, Howell BA, Yang Y, et al. **An analysis of N-acetylcysteine treatment for acetaminophen overdose using a systems model of drug-induced liver injury.** *J Pharmacol Exp Ther.* 2012;342(2):529-40. — prior work explaining the NAC time window with QSP. The direct precedent for analysis [4] of this model. <https://pubmed.ncbi.nlm.nih.gov/22645248/>
100. ★ Woodhead JL, Yang K, Brouwer KLR, et al. **Mechanistic modeling reveals the critical knowledge gaps in bile acid-mediated DILI.** *CPT Pharmacometrics Syst Pharmacol.* 2014;3(3):e123. <https://pubmed.ncbi.nlm.nih.gov/24646538/>
101. ★ Longo DM, Yang Y, Watkins PB, Howell BA, Siler SQ. **Elucidating Differences in the Hepatotoxic Potential of Tolcapone and Entacapone With DILIsym, a Mechanistic Model of Drug-Induced Liver Injury.** *CPT Pharmacometrics Syst Pharmacol.* 2016;5(1):31-9. <https://pubmed.ncbi.nlm.nih.gov/26844013/>
102. ★ Shoda LKM, Woodhead JL, Siler SQ, Watkins PB, Howell BA. **Linking physiology to toxicity using DILIsym®, a mechanistic mathematical model of drug-induced liver injury.** *Biopharm Drug Dispos.* 2014;35(1):33-49. <https://pubmed.ncbi.nlm.nih.gov/24214486/>
103. ★ Baron KT, Gastonguay MR. **Simulation from ODE-based population PK/PD and systems pharmacology models in R with mrgsolve.** *J Pharmacokinet Pharmacodyn.* 2015;42:S84-85. — mrgsolve. <https://github.com/metrumresearchgroup/mrgsolve>
104. Bhattacharya S, Shoda LKM, Zhang Q, et al. **Modeling drug- and chemical-induced hepatotoxicity with systems biology approaches.** *Front Physiol.* 2012;3:462. <https://pubmed.ncbi.nlm.nih.gov/23248599/>

---

## 16. Bistability and threshold behaviour — the theoretical basis for the model's structure

105. ★ Ferrell JE Jr. **Self-perpetuating states in signal transduction: positive feedback, double-negative feedback and bistability.** *Curr Opin Cell Biol.* 2002;14(2):140-8. — the conditions under which a positive feedback loop with saturating gain creates bistability. The interpretive frame for the JNK–Sab loop of the model. <https://pubmed.ncbi.nlm.nih.gov/11891111/>
106. ★ Tyson JJ, Chen KC, Novak B. **Sniffers, buzzers, toggles and blinkers: dynamics of regulatory and signaling pathways in the cell.** *Curr Opin Cell Biol.* 2003;15(2):221-31. <https://pubmed.ncbi.nlm.nih.gov/12648679/>
107. ★ Bagci EZ, Vodovotz Y, Billiar TR, Ermentrout GB, Bahar I. **Bistability in apoptosis: roles of bax, bcl-2, and mitochondrial permeability transition pores.** *Biophys J.* 2006;90(5):1546-59. — a precedent for a bistable model of the cell-death decision that includes the MPT. <https://pubmed.ncbi.nlm.nih.gov/16339882/>

---

### Citation conventions
- PubMed links take the form `https://pubmed.ncbi.nlm.nih.gov/<PMID>/`.
- For the meaning of the ★/▲ marks, see the top of the document.
- Total references: 108 (including a large number of ★ structure-determining references).
- This reference list is for the documentation of a model for educational and research purposes, and does not
  replace clinical practice guidelines.
