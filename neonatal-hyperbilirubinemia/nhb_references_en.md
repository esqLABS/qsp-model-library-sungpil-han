# Neonatal hyperbilirubinaemia QSP model — references
# Neonatal Hyperbilirubinaemia (NHB) QSP Model — References

The evidence for **every structural assumption and parameter** that went into `nhb_qsp_model_en.dot`, `nhb_mrgsolve_model_en.R`, `nhb_shiny_app_en.R`,
and `nhb_reference_check.py` is collected in this document.

**97 references in total** · every PMID was **looked up** through the PubMed E-utilities (esearch + esummary) and the
author, journal, year, and title taken verbatim from it. Not one PMID here was written from memory, and
the bibliographic details were inserted mechanically from the strings the API returned.

The annotation on each entry records "**which parameter or which structural choice** of the model this
paper supports". References not used in the model are not included.

All the computed values behind the quantitative claims come from `nhb_reference_output.txt` (the run
output of an independent python implementation), and A1-A14 in the text refer to that file's section numbers.

---

## 1. Thresholds and clinical guidelines — the model's `thr_pt` / `thr_et` are an analytic approximation to these curves

The AAP 2022 curves are a function of gestational age, time, and neurotoxicity risk factors. The model approximates that structure as `plateau(GA) x (0.42 + 0.58(1-e^{-t/36h}))`, and the risk factors subtract 2 mg/dL from the plateau. **In A1b, albumin of 3.0 g/dL on its own comes out as -2.65 mg/dL — a value reproduced from binding stoichiometry alone, without the guideline ever having been shown to the model.**

1. **Kemper AR** et al. *Clinical Practice Guideline Revision: Management of Hyperbilirubinemia in the Newborn Infant 35 or More Weeks of Gestation.* Pediatrics. 2022.
   [PMID 35927462](https://pubmed.ncbi.nlm.nih.gov/35927462/)  
   → The revised AAP 2022 guideline — the phototherapy and exchange transfusion threshold curves, the list of neurotoxicity risk factors (GA<38 weeks · albumin<3.0 · isoimmune haemolysis · G6PD · sepsis · clinical instability), and permission for home phototherapy. The direct basis for the model's `GA`, `RF`, `thr_pt`, `thr_et`.

2. **Slaughter JL** et al. *Technical Report: Diagnosis and Management of Hyperbilirubinemia in the Newborn Infant 35 or More Weeks of Gestation.* Pediatrics. 2022.
   [PMID 35927519](https://pubmed.ncbi.nlm.nih.gov/35927519/)  
   → The technical report behind the guideline above — the derivation of the thresholds and the evidence grading of the IVIG and exchange transfusion recommendations. The source of the `ESCAL = ET - 2` (escalation-of-care) definition.

3. **Bhutani VK** et al. *Predictive ability of a predischarge hour-specific serum bilirubin for subsequent significant hyperbilirubinemia in healthy term and near-term newborns.* Pediatrics. 1999.
   [PMID 9917432](https://pubmed.ncbi.nlm.nih.gov/9917432/)  
   → The original Bhutani nomogram paper — the pre-discharge hour-specific percentile predicts later severe hyperbilirubinaemia. Why the model treats the time axis as absolute postnatal hours (h).

4. **Bhutani VK** et al. *Noninvasive measurement of total serum bilirubin in a multiracial predischarge newborn population to assess the risk of severe hyperbilirubinemia.* Pediatrics. 2000.
   [PMID 10920173](https://pubmed.ncbi.nlm.nih.gov/10920173/)  
   → The risk-assessment performance of non-invasive total bilirubin measurement in a multi-ethnic cohort — the standing of the model's TcB output (`TCB = 0.92 x C_ex`) as a screening index.

5. **Bhutani VK** et al. *Phototherapy to prevent severe neonatal hyperbilirubinemia in the newborn infant 35 or more weeks of gestation.* Pediatrics. 2011.
   [PMID 21949150](https://pubmed.ncbi.nlm.nih.gov/21949150/)  
   → The phototherapy technical report — quantitative recommendations on irradiance, exposed area, and light source spectrum. The model's `IRRSET`, `FBSASET`, and the definition 'intensive = 30 µW/cm²/nm'.

6. **Wickremasinghe AC** et al. *Neonatal Hyperbilirubinemia.* Pediatr Clin North Am. 2025.
   [PMID 40619190](https://pubmed.ncbi.nlm.nih.gov/40619190/)  
   → The most recent comprehensive review of neonatal hyperbilirubinaemia — the clinical frame for the model's whole structure.

---

## 2. Free bilirubin and albumin binding — the model's first theme

The value that is measured (TSB) and the value that does the damage (Bf) are linked by a saturable binding isotherm. The model solves the **two-site isotherm** — a high-affinity primary site (KA1 ≈ 4.5x10⁷ M⁻¹) and a low-affinity secondary site (KA2 ≈ 1x10⁶ M⁻¹) — by bisection. A single-site isotherm cannot be used, because once the B/A molar ratio exceeds 1 it predicts a physically impossible Bf.

7. **Ahlfors CE** et al. *Unbound (free) bilirubin: improving the paradigm for evaluating neonatal jaundice.* Clin Chem. 2009.
   [PMID 19423734](https://pubmed.ncbi.nlm.nih.gov/19423734/)  
   → The free bilirubin paradigm — why Bf and not TSB is the risk index. Why the model computes `BF` as a separate output and drives BBB influx from Bf alone.

8. **Hulzebos CV** et al. *Bilirubin-albumin binding, bilirubin/albumin ratios, and free bilirubin levels: where do we stand?* Semin Perinatol. 2014.
   [PMID 25304058](https://pubmed.ncbi.nlm.nih.gov/25304058/)  
   → A review setting out where bilirubin–albumin binding, the B/A ratio, and free bilirubin currently stand. The background to the derived result of model A1b, that 'B/A is independent of albumin but not independent of affinity'.

9. **Hulzebos CV** et al. *The bilirubin albumin ratio in the management of hyperbilirubinemia in preterm infants to improve neurodevelopmental outcome: a randomized controlled trial--BARTrial.* PLoS One. 2014.
   [PMID 24927259](https://pubmed.ncbi.nlm.nih.gov/24927259/)  
   → BARTrial — the randomised trial that added the B/A ratio to treatment decisions in preterm infants. Both the basis for using B/A as a clinical index and a demonstration of its limits.

10. **Hegyi T** et al. *Neonatal hyperbilirubinemia and the role of unbound bilirubin.* J Matern Fetal Neonatal Med. 2022.
   [PMID 34957902](https://pubmed.ncbi.nlm.nih.gov/34957902/)  
   → The role of unbound bilirubin — measurement of Bf and its clinical interpretation. The basis for setting the model's Bf injury threshold at 30-35 nM.

11. **Wennberg RP** et al. *Toward understanding kernicterus: a challenge to improve the management of jaundiced newborns.* Pediatrics. 2006.
   [PMID 16452368](https://pubmed.ncbi.nlm.nih.gov/16452368/)  
   → A perspective for understanding kernicterus — risk cannot be defined by TSB alone. The clinical statement of the model's first theme.

12. **Amin SB** et al. *Effect of ibuprofen on bilirubin-albumin binding affinity in premature infants.* J Perinat Med. 2011.
   [PMID 20954849](https://pubmed.ncbi.nlm.nih.gov/20954849/)  
   → Ibuprofen lowers bilirubin–albumin binding affinity in preterm infants — the direct basis for the model's `FDISP` parameter.

13. **Amin SB** et al. *Intravenous lipid and bilirubin-albumin binding variables in premature infants.* Pediatrics. 2009.
   [PMID 19564302](https://pubmed.ncbi.nlm.nih.gov/19564302/)  
   → Intravenous lipid (free fatty acids) changes the binding variables — why a free fatty acid route is included in `FDISP`.

14. **Kleinfeld A** et al. *Displacement of Bilirubin From Albumin by Omegaven and Intralipid.* J Pediatr Surg. 2025.
   [PMID 40946844](https://pubmed.ncbi.nlm.nih.gov/40946844/)  
   → Displacement of bilirubin from albumin by Omegaven/Intralipid — the most recent quantitative basis for `FDISP`.

15. **Martin E** et al. *Ceftriaxone--bilirubin-albumin interactions in the neonate: an in vivo study.* Eur J Pediatr. 1993.
   [PMID 8335024](https://pubmed.ncbi.nlm.nih.gov/8335024/)  
   → An in vivo study of the ceftriaxone–bilirubin–albumin interaction — the `FDISP = 0.70` scenario.

16. **Gulian JM** et al. *Bilirubin displacement by ceftriaxone in neonates: evaluation by determination of 'free' bilirubin and erythrocyte-bound bilirubin.* J Antimicrob Chemother. 1987.
   [PMID 3610909](https://pubmed.ncbi.nlm.nih.gov/3610909/)  
   → Ceftriaxone displacement measured as 'free' bilirubin and erythrocyte-bound bilirubin — direct evidence that displacement raises Bf.

17. **Wadsworth SJ** et al. *In vitro displacement of bilirubin by antibiotics and 2-hydroxybenzoylglycine in newborns.* Antimicrob Agents Chemother. 1988.
   [PMID 3190184](https://pubmed.ncbi.nlm.nih.gov/3190184/)  
   → In vitro quantification of antibiotic-induced bilirubin displacement in neonatal serum — the ranking of displacing potency by drug.

18. **Amin SB** et al. *Unbound Bilirubin and Auditory Neuropathy Spectrum Disorder in Late Preterm and Term Infants with Severe Jaundice.* J Pediatr. 2016.
   [PMID 26952116](https://pubmed.ncbi.nlm.nih.gov/26952116/)  
   → In severely jaundiced term and late preterm infants, unbound bilirubin is associated with auditory neuropathy (ANSD) — the basis for making the model's `BBRTHRA` (ABR threshold) Bf-based.

19. **Ahlfors CE** et al. *Unbound bilirubin predicts abnormal automated auditory brainstem response in a diverse newborn population.* J Perinatol. 2009.
   [PMID 19242487](https://pubmed.ncbi.nlm.nih.gov/19242487/)  
   → Unbound bilirubin predicts abnormal automated ABR (better than TSB does) — the quantitative basis for the Bf → ABR pathway.

---

## 3. Bilirubin production, haem catabolism, and ETCOc — the model's second theme

TSB is the integral of the difference between two fluxes, so on its own it cannot identify the cause. Because HO-1 releases exactly one molecule of CO per molecule of bilirubin, COHb/ETCOc is **a direct read-out of the production flux**. The model implements this as `dCOHB/dt = KCO x PROD/W - KCOOUT x COHB`.

20. **Stevenson DK** et al. *Bilirubin production in healthy term infants as measured by carbon monoxide in breath.* Clin Chem. 1994.
   [PMID 7923775](https://pubmed.ncbi.nlm.nih.gov/7923775/)  
   → Bilirubin production in healthy term infants measured by exhaled CO — the basis for the model's production of 8.5 mg/kg/day (about twice the adult value) and for `BILPERHB = 34 mg/g`, `FEARLY = 0.25`.

21. **Stevenson DK** et al. *Prediction of hyperbilirubinemia in near-term and term infants.* J Perinatol. 2001.
   [PMID 11803421](https://pubmed.ncbi.nlm.nih.gov/11803421/)  
   → Prediction of hyperbilirubinaemia in near-term and term infants — the relation between production indices and the TSB trajectory.

22. **Stevenson DK** et al. *Increased Carbon Monoxide Washout Rates in Newborn Infants.* Neonatology. 2020.
   [PMID 31634890](https://pubmed.ncbi.nlm.nih.gov/31634890/)  
   → The increased rate of CO washout in the newborn — the basis for the model's `KCOOUT` (t½ ≈ 5 h).

23. **Tidmarsh GF** et al. *End-tidal carbon monoxide and hemolysis.* J Perinatol. 2014.
   [PMID 24743136](https://pubmed.ncbi.nlm.nih.gov/24743136/)  
   → End-tidal CO and haemolysis — the basis for ETCOc as the clinical index that separates haemolysis from its absence. The central claim of model A3.

24. **Bhutani VK** et al. *Identification of risk for neonatal haemolysis.* Acta Paediatr. 2018.
   [PMID 29532503](https://pubmed.ncbi.nlm.nih.gov/29532503/)  
   → Identifying haemolytic risk in the newborn — how to distinguish a production lesion from an elimination lesion clinically.

25. **Herschel M** et al. *Isoimmunization is unlikely to be the cause of hemolysis in ABO-incompatible but direct antiglobulin test-negative neonates.* Pediatrics. 2002.
   [PMID 12093957](https://pubmed.ncbi.nlm.nih.gov/12093957/)  
   → In ABO-incompatible but DAT-negative newborns, isoimmunisation is not the cause of haemolysis — why the model puts the haemolytic load in a continuous variable, `ABMAT` (antibody load), rather than in blood group incompatibility itself.

26. **Wong RJ** et al. *In vitro inhibition of heme oxygenase isoenzymes by metalloporphyrins.* J Perinatol. 2011.
   [PMID 21448202](https://pubmed.ncbi.nlm.nih.gov/21448202/)  
   → In vitro inhibition of the HO isoforms by metalloporphyrins — the model's `KISNMP` (competitive inhibition constant) and the site of action of stannsoporfin.

27. **Fujioka K** et al. *Inhibition of heme oxygenase activity using a microparticle formulation of zinc protoporphyrin in an acute hemolytic newborn mouse model.* Pediatr Res. 2016.
   [PMID 26488552](https://pubmed.ncbi.nlm.nih.gov/26488552/)  
   → Inhibition of HO activity by a microparticle formulation of zinc protoporphyrin in acutely haemolysing newborn mice — causal proof that blocking production lowers bilirubin.

---

## 4. UGT1A1 maturation · genotype · hepatic uptake — the elimination flux

The model writes this as `UGT_target = GENO x ont(t) x (1 + E_pb C/(EC50+C)) + TGX`. In that multiplicative structure, **if GENO = 0 (CN-I) the induction term vanishes with it**, and that produces — with no separate rule — the reason phenobarbital works in CN-II and not in CN-I.

28. **Nie YL** et al. *Histone Modifications Regulate the Developmental Expression of Human Hepatic UDP-Glucuronosyltransferase 1A1.* Drug Metab Dispos. 2017.
   [PMID 29025858](https://pubmed.ncbi.nlm.nih.gov/29025858/)  
   → Developmental expression of UGT1A1 in human liver is regulated by histone modification — the molecular basis for the postnatal maturation time constant `TAUUGT`.

29. **Nie YL** et al. *Hepatic expression of transcription factors affecting developmental regulation of UGT1A1 in the Han Chinese population.* Eur J Clin Pharmacol. 2017.
   [PMID 27704169](https://pubmed.ncbi.nlm.nih.gov/27704169/)  
   → Hepatic expression of the transcription factors regulating UGT1A1 development in a Han Chinese population — population variation in the maturation curve.

30. **Travan L** et al. *Severe neonatal hyperbilirubinemia and UGT1A1 promoter polymorphism.* J Pediatr. 2014.
   [PMID 24726540](https://pubmed.ncbi.nlm.nih.gov/24726540/)  
   → Severe neonatal hyperbilirubinaemia and UGT1A1 promoter polymorphism — the clinical basis for `GENO` (Gilbert 0.35).

31. **Žaja O** et al. *Correlation of UGT1A1 TATA-box polymorphism and jaundice in breastfed newborns-early presentation of Gilbert's syndrome.* J Matern Fetal Neonatal Med. 2014.
   [PMID 23981182](https://pubmed.ncbi.nlm.nih.gov/23981182/)  
   → The correlation between UGT1A1 TATA-box polymorphism and breastfeeding jaundice in the newborn — the genotype x feeding interaction of model A9.

32. **Yamamoto A** et al. *Gly71Arg mutation of the bilirubin UDP-glucuronosyltransferase 1A1 gene is associated with neonatal hyperbilirubinemia in the Japanese population.* Kobe J Med Sci. 2002.
   [PMID 12502904](https://pubmed.ncbi.nlm.nih.gov/12502904/)  
   → The association of Gly71Arg (UGT1A1*6) with neonatal hyperbilirubinaemia in a Japanese population — `GENO = 0.40` (the East Asian variant).

33. **Prachukthum S** et al. *Association between UGT 1A1 Gly71Arg (G71R) polymorphism and neonatal hyperbilirubinemia.* J Med Assoc Thai. 2012.
   [PMID 23964438](https://pubmed.ncbi.nlm.nih.gov/23964438/)  
   → UGT1A1 G71R and neonatal hyperbilirubinaemia in Thai infants — replication of the East and Southeast Asian variant.

34. **Talebi H** et al. *Association of SLCO1B1 genetic variants with neonatal hyperbilirubinemia: a consolidated analysis of 36 studies.* BMC Pediatr. 2025.
   [PMID 40155882](https://pubmed.ncbi.nlm.nih.gov/40155882/)  
   → SLCO1B1 variants and neonatal hyperbilirubinaemia pooled across 36 studies — the model's `FOATP` (hepatic uptake genotype coefficient).

35. **Fan J** et al. *Associations between UGT1A1, SLCO1B1, SLCO1B3, BLVRA and HMOX1 polymorphisms and susceptibility to neonatal severe hyperbilirubinemia in Chinese Han population.* BMC Pediatr. 2024.
   [PMID 38279097](https://pubmed.ncbi.nlm.nih.gov/38279097/)  
   → UGT1A1·SLCO1B1·SLCO1B3·BLVRA·HMOX1 polymorphisms and severe hyperbilirubinaemia — the basis for the model placing a genetic coefficient at each of the uptake, conjugation, and production steps.

36. **Wang WF** et al. *Variants in UGT1A1 and SLCO1B1 increase the risk of neonatal hyperbilirubinemia: a case-control study in subtropical China.* Front Pediatr. 2026.
   [PMID 42422451](https://pubmed.ncbi.nlm.nih.gov/42422451/)  
   → The increased risk from UGT1A1 and SLCO1B1 variants in a subtropical Chinese cohort — the basis for uptake and conjugation being independent risk axes.

---

## 5. Crigler-Najjar and enzyme replacement — the model's fourth theme

Because lumirubin requires no conjugation, phototherapy is the only route that bypasses UGT1A1 completely. The model computes that about 10 % of adult UGT1A1 activity is enough to turn CN-I into a CN-II phenotype, with TSB staying at ~9 mg/dL without lamps (A11c).

37. **D'Antiga L** et al. *Gene Therapy in Patients with the Crigler-Najjar Syndrome.* N Engl J Med. 2023.
   [PMID 37585628](https://pubmed.ncbi.nlm.nih.gov/37585628/)  
   → GNT0003, the first human gene therapy — after AAV8-hUGT1A1 in CN-I patients, bilirubin was controlled without phototherapy. The direct basis for the model's `TGXMAX = 0.10` and `KTGON` (t½ 10 days).

38. **Collaud F** et al. *Preclinical Development of an AAV8-hUGT1A1 Vector for the Treatment of Crigler-Najjar Syndrome.* Mol Ther Methods Clin Dev. 2019.
   [PMID 30705921](https://pubmed.ncbi.nlm.nih.gov/30705921/)  
   → Preclinical development of the AAV8-hUGT1A1 vector — expression level and dose-response.

39. **Aronson SJ** et al. *Prevalence and Relevance of Pre-Existing Anti-Adeno-Associated Virus Immunity in the Context of Gene Therapy for Crigler-Najjar Syndrome.* Hum Gene Ther. 2019.
   [PMID 31502485](https://pubmed.ncbi.nlm.nih.gov/31502485/)  
   → The prevalence and significance of pre-existing anti-AAV immunity in CN gene therapy — cited to make explicit a failure mode the model does not cover (neutralising antibodies).

40. **Agati G** et al. *Bilirubin photoisomerization products in serum and urine from a Crigler-Najjar type I patient treated by phototherapy.* J Photochem Photobiol B. 1998.
   [PMID 10093917](https://pubmed.ncbi.nlm.nih.gov/10093917/)  
   → Serum and urine photoisomers in CN-I patients receiving phototherapy — **direct evidence that the photochemical products really are excreted in humans, and the decisive observation that phototherapy works independently of conjugation.** The experimental basis for the model's third and fourth themes.

41. **Jesus Sá C** et al. *Between Crigler-Najjar Syndrome Type II and Gilbert Syndrome: Expanding the Spectrum of Uridine Diphosphate Glucuronosyltransferase 1A1 (UGT1A1)-Related Hyperbilirubinemia.* Cureus. 2026.
   [PMID 41658774](https://pubmed.ncbi.nlm.nih.gov/41658774/)  
   → The spectrum between CN-II and Gilbert — the basis for `GENO` as a continuous activity scale (CN-I 0 → CN-II 0.05 → Gilbert 0.35 → wild type 1.0).

---

## 6. Photochemistry and the phototherapy dose-response — the model's third theme

Phototherapy is two reactions. A fast, **reversible** (4Z,15Z) ⇄ (4Z,15E) configurational isomerisation reaches a photostationary state (20.6 % of measured TSB in the model), and an **irreversible** cyclisation of the E-isomer forms lumirubin, which does the actual excreting. Both photoisomers are measured as 'bilirubin' by routine assays, so 26 % of the reported TSB is photoisomer; and when the lamps are switched off, TSB first **falls** (lumirubin washout, t½ 1 h) and then rises again through thermal back-isomerisation.

42. **Tan KL** *The nature of the dose-response relationship of phototherapy for neonatal hyperbilirubinemia.* J Pediatr. 1977.
   [PMID 839340](https://pubmed.ncbi.nlm.nih.gov/839340/)  
   → The nature of the phototherapy dose-response relationship — diminishing returns as irradiance is raised. The classic basis for the model's `I50 = 45 µW/cm²/nm` (optical saturation).

43. **Tan KL** *The pattern of bilirubin response to phototherapy for neonatal hyperbilirubinaemia.* Pediatr Res. 1982.
   [PMID 7110789](https://pubmed.ncbi.nlm.nih.gov/7110789/)  
   → Patterns of the bilirubin response to phototherapy — the size of the 24-hour fall (model: 39 % under intensive phototherapy).

44. **Pratesi R** et al. *Phototherapy for neonatal hyperbilirubinemia.* Photodermatol. 1989.
   [PMID 2700092](https://pubmed.ncbi.nlm.nih.gov/2700092/)  
   → Phototherapy for neonatal hyperbilirubinaemia — a survey of light source spectra, penetration depth, and photochemical routes. The basis for the model decomposing photon delivery as `irradiance x exposed area x optical accessibility`.

45. **Borden AR** et al. *Variation in the Phototherapy Practices and Irradiance of Devices in a Major Metropolitan Area.* Neonatology. 2018.
   [PMID 29393277](https://pubmed.ncbi.nlm.nih.gov/29393277/)  
   → Phototherapy practice and the variation in actual irradiance across metropolitan-area equipment — the point that the prescribed irradiance and the delivered irradiance differ. Why the model states that its `IRRSET` is the delivered irradiance.

46. **Chang PW** et al. *A Clinical Prediction Rule for Rebound Hyperbilirubinemia Following Inpatient Phototherapy.* Pediatrics. 2017.
   [PMID 28196932](https://pubmed.ncbi.nlm.nih.gov/28196932/)  
   → A clinical rule for predicting rebound after phototherapy is stopped — the comparator for the reading that the model's rebound (A5: +2.2 mg/dL) is in part photochemical bookkeeping (E→Z back-conversion).

47. **Slusher TM** et al. *A Randomized Trial of Phototherapy with Filtered Sunlight in African Neonates.* N Engl J Med. 2015.
   [PMID 26376136](https://pubmed.ncbi.nlm.nih.gov/26376136/)  
   → A randomised trial of filtered sunlight phototherapy — evidence from a resource-limited setting that what matters is photon delivery and not the power source.

48. **Slusher TM** et al. *Filtered sunlight versus intensive electric powered phototherapy in moderate-to-severe neonatal hyperbilirubinaemia: a randomised controlled non-inferiority trial.* Lancet Glob Health. 2018.
   [PMID 30170894](https://pubmed.ncbi.nlm.nih.gov/30170894/)  
   → A non-inferiority trial of filtered sunlight vs intensive electric phototherapy (moderate-severe) — as above.

49. **Magee S** et al. *Multi-directional phototherapy vs. unidirectional phototherapy from below for severe neonatal jaundice: a randomized pilot trial in home phototherapy.* Eur J Pediatr. 2026.
   [PMID 41975040](https://pubmed.ncbi.nlm.nih.gov/41975040/)  
   → A randomised trial of multidirectional vs single-direction downward phototherapy — **the clinical trial corresponding to the conclusion of model A4, that exposed body surface area matters more than irradiance.**

50. **Cordero N** et al. *Transcutaneous bilirubin measurements in preterm infants: the impact of race, age, and phototherapy.* J Perinatol. 2026.
   [PMID 41490938](https://pubmed.ncbi.nlm.nih.gov/41490938/)  
   → The effect of ethnicity, age, and phototherapy on TcB in preterm infants — why the model computes TcB from the skin compartment (`BEX`) and lets it dissociate from serum under the lamps.

---

## 7. The enterohepatic circulation — the newborn's second liver, but running the other way

In the newborn, conjugated bilirubin is not an exit but a **recirculation**. β-glucuronidase deconjugates it, the jejunum reabsorbs it, and the portal vein returns it to the liver. The effective clearance is therefore `CL_conj x (1 - f_reabsorb)`, and feeding support carries the same units as a drug effect.

51. **Kreamer BL** et al. *A novel inhibitor of beta-glucuronidase: L-aspartic acid.* Pediatr Res. 2001.
   [PMID 11568288](https://pubmed.ncbi.nlm.nih.gov/11568288/)  
   → A β-glucuronidase inhibitor (L-aspartic acid) — proof that the deconjugation step can be a pharmacological target. The scope for intervening on the model's `KBGLUC` term.

52. **Alonso EM** et al. *Enterohepatic circulation of nonconjugated bilirubin in rats fed with human milk.* J Pediatr. 1991.
   [PMID 1999786](https://pubmed.ncbi.nlm.nih.gov/1999786/)  
   → Enterohepatic circulation of unconjugated bilirubin in breast-milk-fed rats — the basis for the `BGABM` term (increased β-glucuronidase from breast milk).

53. **Yau KI** et al. *Factors affecting the severity of neonatal jaundice of unknown etiology: the role of enterohepatic circulation.* Zhonghua Min Guo Xiao Er Ke Yi Xue Hui Za Zhi. 1992.
   [PMID 1626448](https://pubmed.ncbi.nlm.nih.gov/1626448/)  
   → The enterohepatic circulation contributes to the severity of neonatal jaundice of unknown cause — the basis for a recirculated fraction of ~50 %.

54. **Vítek L** et al. *The impact of intestinal microflora on serum bilirubin levels.* J Hepatol. 2005.
   [PMID 15664250](https://pubmed.ncbi.nlm.nih.gov/15664250/)  
   → The effect of the gut microbiota on serum bilirubin levels — the point that a germ-free gut performs no urobilinoid conversion. The model's `TAUBGA` (flora maturation 21 days).

55. **Vítek L** et al. *Gut microbiota and bilirubin metabolism: unveiling new pathways in health and disease.* Trends Mol Med. 2025.
   [PMID 39757046](https://pubmed.ncbi.nlm.nih.gov/39757046/)  
   → The most recent review of the gut microbiota and bilirubin metabolism — a modern account of the pathway above.

56. **Vítek L** et al. *Identification of bilirubin reduction products formed by Clostridium perfringens isolated from human neonatal fecal flora.* J Chromatogr B Analyt Technol Biomed Life Sci. 2006.
   [PMID 16504607](https://pubmed.ncbi.nlm.nih.gov/16504607/)  
   → Bilirubin reduction products formed by Clostridium perfringens in neonatal faecal flora — the mechanism by which establishment of the flora breaks the recirculation.

57. **Abdel-Aziz Ali SM** et al. *Efficacy of oral agar in management of indirect hyperbilirubinemia in full-term neonates.* J Matern Fetal Neonatal Med. 2022.
   [PMID 32192396](https://pubmed.ncbi.nlm.nih.gov/32192396/)  
   → The efficacy of oral agar in indirect hyperbilirubinaemia of term infants — the clinical basis for the model's `GBIND`/`KB50` (enteral binder) terms.

58. **Lazarus G** et al. *Role of ursodeoxycholic acid in neonatal indirect hyperbilirubinemia: a systematic review and meta-analysis of randomized controlled trials.* Ital J Pediatr. 2022.
   [PMID 36253867](https://pubmed.ncbi.nlm.nih.gov/36253867/)  
   → The role of UDCA in neonatal indirect hyperbilirubinaemia: a meta-analysis of randomised trials — `EMAXU`, `EC50UDCA`.

59. **Tabrizi M** et al. *Ursodeoxycholic Acid in the Management of Prolonged Neonatal Hyperbilirubinemia: A Randomized Controlled Clinical Trial.* Iran J Med Sci. 2026.
   [PMID 42100234](https://pubmed.ncbi.nlm.nih.gov/42100234/)  
   → A randomised trial of UDCA in prolonged neonatal hyperbilirubinaemia — as above, in prolonged jaundice.

60. **de Oliveira HM** et al. *Zinc sulfate on neonatal hyperbilirubinemia: an updated systematic review and meta-analysis.* Eur J Pediatr. 2024.
   [PMID 39671002](https://pubmed.ncbi.nlm.nih.gov/39671002/)  
   → An updated systematic review and meta-analysis of zinc sulphate — blocking reabsorption by intraluminal precipitation (the `GBIND` family in the model).

61. **Gholitabar M** et al. *Clofibrate in combination with phototherapy for unconjugated neonatal hyperbilirubinaemia.* Cochrane Database Syst Rev. 2012.
   [PMID 23235669](https://pubmed.ncbi.nlm.nih.gov/23235669/)  
   → The Cochrane review of clofibrate plus phototherapy — UGT1A1 induction through PPARα (an alternative agonist for the model's `CAR` route).

62. **Fei Q** et al. *Phototherapy Alone or Combined with Adjuvant Drugs for Neonatal Hyperbilirubinemia: A Systematic Review and Network Meta-Analysis.* Children (Basel). 2026.
   [PMID 42073150](https://pubmed.ncbi.nlm.nih.gov/42073150/)  
   → A network meta-analysis of phototherapy alone vs phototherapy with adjuvant drugs — the most recent quantitative standard against which to compare the model's ranking of adjuvants (A10).

63. **Tavares LC** et al. *Oral adjuvants supplemented to phototherapy as a therapeutic strategy for neonatal hyperbilirubinemia: a systematic review.* J Perinat Med. 2026.
   [PMID 41811691](https://pubmed.ncbi.nlm.nih.gov/41811691/)  
   → A systematic review of oral drugs adjuvant to phototherapy — as above.

---

## 8. Exchange transfusion and IVIG — acting on load, binding capacity, and production at once

In the model an exchange transfusion is not an instantaneous event but a real procedure spread over 2-4 hours, and the amount removed **comes out** of the integral `∫ Q_ET x C_p dt`. Because donor plasma resets albumin and removes maternal antibody, the fall in Bf (82 %) is larger than the fall in TSB (40 %).

64. **Zwiers C** et al. *Immunoglobulin for alloimmune hemolytic disease in neonates.* Cochrane Database Syst Rev. 2018.
   [PMID 29551014](https://pubmed.ncbi.nlm.nih.gov/29551014/)  
   → The Cochrane review of immunoglobulin in neonatal isoimmune haemolytic disease — the basis for `IMAXIVIG = 0.65` and `IC50IVIG`, and an explicit statement of the uncertainty in the effect size.

65. **Huang M** et al. *Efficacy and safety of different doses of intravenous immunoglobulin combined with phototherapy in neonatal hemolytic disease: a systematic review and meta-analysis.* Eur J Pediatr. 2026.
   [PMID 42458120](https://pubmed.ncbi.nlm.nih.gov/42458120/)  
   → A meta-analysis of efficacy and safety by IVIG dose — the basis for choosing the 1 g/kg dose.

66. **van Klink JM** et al. *Immunoglobulins in Neonates with Rhesus Hemolytic Disease of the Fetus and Newborn: Long-Term Outcome in a Randomized Trial.* Fetal Diagn Ther. 2016.
   [PMID 26159803](https://pubmed.ncbi.nlm.nih.gov/26159803/)  
   → Long-term outcome of immunoglobulin in newborns with Rh haemolytic disease (a randomised trial) — contrary evidence on whether IVIG reduces the need for exchange transfusion. Why the model caps the IVIG effect at an Emax of 0.65.

67. **Legler TJ** *RhIg for the prevention Rh immunization and IVIg for the treatment of affected neonates.* Transfus Apher Sci. 2020.
   [PMID 33004277](https://pubmed.ncbi.nlm.nih.gov/33004277/)  
   → Prophylactic anti-D (RhIg) and therapeutic IVIg — the basis for placing RhIG as 'primary prevention' on the map.

68. **Dash N** et al. *Pre exchange Albumin Administration in Neonates with Hyperbilirubinemia: A Randomized Controlled Trial.* Indian Pediatr. 2015.
   [PMID 26519710](https://pubmed.ncbi.nlm.nih.gov/26519710/)  
   → A randomised trial of albumin before exchange transfusion — **clinical proof that raising binding capacity is a treatment axis separate from lowering the load.** The model's `ALBINF`/`ALBDON` terms.

69. **Shahian M** et al. *Effect of albumin administration prior to exchange transfusion in term neonates with hyperbilirubinemia--a randomized controlled trial.* Indian Pediatr. 2010.
   [PMID 19578230](https://pubmed.ncbi.nlm.nih.gov/19578230/)  
   → A randomised trial of albumin before exchange transfusion (term infants) — as above.

70. **Hemmati F** et al. *Exchange Transfusion Trends and Risk Factors for Extreme Neonatal Hyperbilirubinemia over 10 Years in Shiraz, Iran.* Iran J Med Sci. 2024.
   [PMID 38952637](https://pubmed.ncbi.nlm.nih.gov/38952637/)  
   → Ten-year trends in exchange transfusion for extreme hyperbilirubinaemia and its risk factors — the characteristics of the population that actually receives exchange transfusion.

71. **Sarı EE** et al. *The '72-Hour Divergence' in severe hyperbilirubinemia: sepsis is associated with altered renal and hematologic adaptation in a 10-year exchange transfusion cohort.* BMC Pediatr. 2026.
   [PMID 41776517](https://pubmed.ncbi.nlm.nih.gov/41776517/)  
   → The '72-hour divide' in severe hyperbilirubinaemia — in a ten-year exchange transfusion cohort, sepsis alters the renal and haematological adaptation. The clinical counterpart of the model's `FBBB` (increased barrier permeability from sepsis) and `FACID`.

---

## 9. Haem oxygenase inhibition — production blockade as a separate treatment axis

Stannsoporfin is a competitive inhibitor of HO-1 and acts in the model **on the production term only**, as `F_SNMP = 1/(1 + C/KISNMP)`. That is why it is meaningful only in a production lesion, the kind ETCOc identifies, and meaningless in an elimination lesion (CN, Gilbert).

72. **Rosenfeld WN** et al. *Stannsoporfin with phototherapy to treat hyperbilirubinemia in newborn hemolytic disease.* J Perinatol. 2022.
   [PMID 34635771](https://pubmed.ncbi.nlm.nih.gov/34635771/)  
   → A clinical trial of stannsoporfin plus phototherapy in neonatal haemolytic disease — the single 4.5 mg/kg IM dose and the effect size (the model's `KASNMP`, `VSNMP`, `KESNMP`, `KISNMP`).

73. **Bhutani VK** et al. *Clinical trial of tin mesoporphyrin to prevent neonatal hyperbilirubinemia.* J Perinatol. 2016.
   [PMID 26938918](https://pubmed.ncbi.nlm.nih.gov/26938918/)  
   → A clinical trial of tin mesoporphyrin for the prevention of neonatal hyperbilirubinaemia — the basis for a 40-75 % production blockade.

74. **Poudel P** et al. *Efficacy and Safety Concerns with Sn-Mesoporphyrin as an Adjunct Therapy in Neonatal Hyperbilirubinemia: A Literature Review.* Int J Pediatr. 2022.
   [PMID 35898803](https://pubmed.ncbi.nlm.nih.gov/35898803/)  
   → A literature review of the efficacy and safety of Sn-mesoporphyrin as adjuvant therapy — making its unapproved status and the safety issues explicit.

---

## 10. Neurotoxicity — what crosses the barrier and what leaves damage behind

The model lets **only free bilirubin** cross the barrier, as `dBBR/dt = KINBBB x FBBB x Bf - KOUTBBB x BBR`, and accumulates damage in proportion to the excess above threshold. The ABR deficit is separated out as reversible and the accumulated damage as irreversible.

75. **Watchko JF** et al. *Bilirubin-induced neurologic damage--mechanisms and management approaches.* N Engl J Med. 2013.
   [PMID 24256380](https://pubmed.ncbi.nlm.nih.gov/24256380/)  
   → Bilirubin-induced neurological damage — mechanism and management (an authoritative review). The structural source for the whole of the model's seventh cluster.

76. **Ostrow JD** et al. *Molecular basis of bilirubin-induced neurotoxicity.* Trends Mol Med. 2004.
   [PMID 15102359](https://pubmed.ncbi.nlm.nih.gov/15102359/)  
   → The molecular basis of bilirubin neurotoxicity — mitochondrial uncoupling, oxidative stress, and cell death. The mechanism of the model's `KINJ` pathway.

77. **Watchko JF** et al. *Brain bilirubin content is increased in P-glycoprotein-deficient transgenic null mutant mice.* Pediatr Res. 1998.
   [PMID 9803459](https://pubmed.ncbi.nlm.nih.gov/9803459/)  
   → Increased brain bilirubin content in P-glycoprotein-deficient mice — **direct evidence that P-gp efflux belongs in the model's `KOUTBBB`.**

78. **Hansen TW** et al. *Rates of bilirubin clearance from rat brain regions.* Biol Neonate. 1995.
   [PMID 8534773](https://pubmed.ncbi.nlm.nih.gov/8534773/)  
   → Region-specific rates of bilirubin clearance from rat brain — the quantitative basis for `KOUTBBB` (t½ ≈ 7 h), and for brain accumulation being potentially reversible.

79. **Amin SB** et al. *Bilirubin and serial auditory brainstem responses in premature infants.* Pediatrics. 2001.
   [PMID 11335741](https://pubmed.ncbi.nlm.nih.gov/11335741/)  
   → Bilirubin and serial ABR in preterm infants — ABR latency lengthens as Bf rises and recovers after treatment. The basis for the model's **reversible** `ABRD` state.

80. **Amin SB** *Clinical assessment of bilirubin-induced neurotoxicity in premature infants.* Semin Perinatol. 2004.
   [PMID 15686265](https://pubmed.ncbi.nlm.nih.gov/15686265/)  
   → Clinical assessment of bilirubin neurotoxicity in preterm infants — BIND-type scores and their dependence on gestational age.

81. **Alexandra Brito M** et al. *Bilirubin toxicity to human erythrocytes: a review.* Clin Chim Acta. 2006.
   [PMID 16887110](https://pubmed.ncbi.nlm.nih.gov/16887110/)  
   → A review of bilirubin toxicity to erythrocytes — the mechanism shared at the membrane and mitochondrial level.

82. **Lu Y** et al. *Bilirubin Oxidation End Products (BOXes) Induce Neuronal Oxidative Stress Involving the Nrf2 Pathway.* Oxid Med Cell Longev. 2021.
   [PMID 34373769](https://pubmed.ncbi.nlm.nih.gov/34373769/)  
   → Bilirubin oxidation end products (BOXes) induce neural oxidative stress through the Nrf2 pathway — cited as a warning that the model's `KALT` (the non-UGT oxidative degradation route) may not be a harmless exit.

83. **Seidel RA** et al. *Impact of higher-order heme degradation products on hepatic function and hemodynamics.* J Hepatol. 2017.
   [PMID 28412296](https://pubmed.ncbi.nlm.nih.gov/28412296/)  
   → The effect of higher-order haem degradation products on liver function and haemodynamics — as above.

84. **Wisnowski JL** et al. *Magnetic resonance imaging of bilirubin encephalopathy: current limitations and future promise.* Semin Perinatol. 2014.
   [PMID 25267277](https://pubmed.ncbi.nlm.nih.gov/25267277/)  
   → MRI in bilirubin encephalopathy — the limits and prospects of globus pallidus signal change.

85. **Gburek-Augustat J** et al. *Acute and Chronic Kernicterus: MR Imaging Evolution of Globus Pallidus Signal Change during Childhood.* AJNR Am J Neuroradiol. 2023.
   [PMID 37620154](https://pubmed.ncbi.nlm.nih.gov/37620154/)  
   → The evolution through childhood of MRI globus pallidus signal change in acute and chronic kernicterus — the imaging separation of the acute reversible component from the chronic irreversible one.

86. **Yi M** et al. *Globus pallidus/putamen T(1)WI signal intensity ratio in grading and predicting prognosis of neonatal acute bilirubin encephalopathy.* Front Pediatr. 2023.
   [PMID 37842026](https://pubmed.ncbi.nlm.nih.gov/37842026/)  
   → Grading of the globus pallidus/striatum T1WI signal ratio and prognostic prediction — the clinical counterpart of the `INJ` index.

87. **Gelineau-Morel R** et al. *Predictive and diagnostic measures for kernicterus spectrum disorder: a prospective cohort study.* Pediatr Res. 2024.
   [PMID 37689774](https://pubmed.ncbi.nlm.nih.gov/37689774/)  
   → Predictive and diagnostic markers of kernicterus spectrum disorder (a prospective cohort) — the comparator for the `KERN` logistic output.

88. **Johnson L** et al. *Clinical report from the pilot USA Kernicterus Registry (1992 to 2004).* J Perinatol. 2009.
   [PMID 19177057](https://pubmed.ncbi.nlm.nih.gov/19177057/)  
   → The clinical report of the US kernicterus registry (1992-2004) — the TSB distribution and the risk factor combinations of the cases that actually occurred.

89. **Bhutani VK** et al. *Synopsis report from the pilot USA Kernicterus Registry.* J Perinatol. 2009.
   [PMID 19177058](https://pubmed.ncbi.nlm.nih.gov/19177058/)  
   → The overview report of the US kernicterus registry — as above.

---

## 11. G6PD deficiency — the commonest cause of kernicterus worldwide

The stoichiometric point of model A8: destruction of 1 g/dL of haemoglobin yields 29 mg/kg of bilirubin, that is **a TSB of 11.6 mg/dL** in a distribution space of 2.5 dL/kg. A fall of 3 g/dL in Hb therefore delivers about 35 mg/dL of pigment, and that is the quantitative reason bilirubin reaches exchange transfusion levels while the haemoglobin still looks nearly normal.

90. **Howes RE** et al. *Spatial distribution of G6PD deficiency variants across malaria-endemic regions.* Malar J. 2013.
   [PMID 24228846](https://pubmed.ncbi.nlm.nih.gov/24228846/)  
   → The distribution of G6PD deficiency variants in malaria-endemic regions — the background to the `G6PD` flag standing for hundreds of millions of people worldwide.

91. **Bancone G** et al. *Contribution of genetic factors to high rates of neonatal hyperbilirubinaemia on the Thailand-Myanmar border.* PLOS Glob Public Health. 2022.
   [PMID 36962413](https://pubmed.ncbi.nlm.nih.gov/36962413/)  
   → The contribution of genetic factors to the high rate of neonatal hyperbilirubinaemia on the Thailand–Myanmar border — the combined effect of G6PD and UGT1A1 variants.

92. **Kaplan M** et al. *Grand Rounds Hyperbilirubinemia following Phototherapy in Glucose-6-Phosphate Dehydrogenase-Deficient Neonates: Not Out of the Woods.* J Pediatr. 2023.
   [PMID 37169338](https://pubmed.ncbi.nlm.nih.gov/37169338/)  
   → Hyperbilirubinaemia after phototherapy in G6PD-deficient newborns — the clinical observation that the risk remains once phototherapy has been finished. The reading of the model's rebound (A5) and closed-loop phototherapy (A13).

93. **Somanath S** et al. *Kernicterus in a late preterm infant with overlooked G6PD deficiency.* BMJ Case Rep. 2025.
   [PMID 40876910](https://pubmed.ncbi.nlm.nih.gov/40876910/)  
   → Kernicterus in a late preterm infant with overlooked G6PD deficiency — the clinical counterpart of the model's S6·S7.

94. **Zahedpasha Y** et al. *Relation between Neonatal Icter and Gilbert Syndrome in Gloucose-6-Phosphate Dehydrogenase Deficient Subjects.* J Clin Diagn Res. 2014.
   [PMID 24783083](https://pubmed.ncbi.nlm.nih.gov/24783083/)  
   → The relation between neonatal jaundice and Gilbert's syndrome in G6PD-deficient individuals — **a multiplicative interaction between a production lesion and an elimination lesion** (the genotype x load structure of model A9).

95. **Zakerihamidi M** et al. *Comparison of severity and prognosis of jaundice due to Rh incompatibility and G6PD deficiency.* Transfus Apher Sci. 2023.
   [PMID 37164807](https://pubmed.ncbi.nlm.nih.gov/37164807/)  
   → A comparison of severity and prognosis between jaundice from Rh incompatibility and from G6PD deficiency — why the same haemolytic load produces different trajectories (model A8: pulse vs sustained).

96. **Xu JX** et al. *Etiology analysis for term newborns with severe hyperbilirubinemia in eastern Guangdong of China.* World J Clin Cases. 2023.
   [PMID 37123300](https://pubmed.ncbi.nlm.nih.gov/37123300/)  
   → An analysis of the causes of severe hyperbilirubinaemia in term infants in Guangdong, eastern China — the actual distribution of the contribution of each cause.

---

## 12. Body surface area and growth — the quantitative basis of the fourth theme

Every photoreaction rate in the model is `photon delivery x concentration within the irradiated layer`, and never `x amount`. Photon delivery therefore scales with body surface area, and the photo-clearance per kg falls as `BSA/W ~ W^-0.30`. **A11b, however, shows that this term is almost exactly offset by the fall in bilirubin production per kg (7.7 → 3.4 mg/kg/day)** — that is, the usual explanation that 'the lamps stop working because the surface-to-mass ratio shrinks' cannot on its own account for the clinical course.

97. **Haycock GB** et al. *Geometric method for measuring body surface area: a height-weight formula validated in infants, children, and adults.* J Pediatr. 1978.
   [PMID 650346](https://pubmed.ncbi.nlm.nih.gov/650346/)  
   → The Haycock body surface area formula (validated in newborns, children, and adults) — the source of the model's `bsa_cm2()`. BSA = 0.024265 x W^0.5378 x H^0.3964.

---

## What the model deliberately does not cover

The literature exists for the items below, but they were **not put into** the model. They are the boundary within which the results should be read.

- **A deep, slow tissue compartment.** Because the extravascular space is treated as a single homogeneous compartment, a double-volume
  exchange transfusion removes 40-45 % of the total body load (the classic teaching figure is ~25 %). For the same reason
  the phototherapy limit for CN-I is optimistic (4.9-8.4 vs the reported 15-25 mg/dL). The two biases
  have **the same cause** and cannot both be fixed at once with `PSXKG` alone.
- Gene therapy failure from **anti-AAV neutralising antibodies** (31502485 above).
- The cholestatic mechanism of **bronze baby syndrome** is on the map but exists in the ODEs only as the `FCHOL` coefficient.
- The behavioural consequences of **mother–infant separation and interrupted feeding**, and the retinal and DNA oxidative damage of phototherapy, are on the map only.
- **TcB bias by ethnicity and skin pigment** is not reflected in the outputs (41490938 above).
- **The AAP 2022 thresholds are an analytic approximation**, not tabulated values. They cannot be used for clinical judgement.

---

## ⚠️ Disclaimer

This is a qualitative and semi-quantitative QSP model for educational and research purposes. It was built from the
published literature but has not been independently verified or certified, and **must not be used directly for real clinical**
**decision-making, prescribing, or regulatory submission.**
