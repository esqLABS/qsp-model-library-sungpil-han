# Alkaptonuria (AKU) — references
### References for the AKU QSP model

Every citation has had its authors, year, journal and title **verified mechanically** through the
PubMed E-utilities (`esummary`). Not a single PMID written from memory has been left in. The text states
which term of the model (parameter · structure · validation target) each citation was used for, and
**the points at which the model disagrees with the literature are written out under that citation too.**

**131 papers** in total. The distinction between anchors (values used in calibration) and held-out
values (used only for validation) corresponds to the validation table in [README.md](../../../alkaptonuria/README.md).

---

## 1. Definition and comprehensive reviews

1. Bernardini G et al. (2024). *Alkaptonuria.* Nat Rev Dis Primers. [PMID 38453957](https://pubmed.ncbi.nlm.nih.gov/38453957/)
   > The definitive review at present. Primary basis for this model's pathophysiological skeleton (HGD deficiency → HGA accumulation → BQA → collagen polymerisation → ochronosis) and for the distribution of involvement across organs.

2. Davison AS et al. (2023). *Alkaptonuria - Past, present and future.* Adv Clin Chem. [PMID 37268334](https://pubmed.ncbi.nlm.nih.gov/37268334/)
   > A survey of the history, present and future of AKU research from a clinical-chemistry standpoint. Basis for the list of measurable biomarkers.

3. Zatkova A et al. (2020). *Alkaptonuria: Current Perspectives.* Appl Clin Genet. [PMID 32158253](https://pubmed.ncbi.nlm.nih.gov/32158253/)

4. Ranganath LR et al. (2013). *Recent advances in management of alkaptonuria (invited review; best practice article).* J Clin Pathol. [PMID 23486607](https://pubmed.ncbi.nlm.nih.gov/23486607/)

5. Gao Y (2026). *Alkaptonuria.* J Clin Rheumatol. [PMID 41451484](https://pubmed.ncbi.nlm.nih.gov/41451484/)

6. Efridi W et al. (2026). *Ochronosis.* StatPearls. [PMID 32809369](https://pubmed.ncbi.nlm.nih.gov/32809369/)


## 2. HGD genetics and population genetics

7. Zatkova A (2011). *An update on molecular genetics of Alkaptonuria (AKU).* J Inherit Metab Dis. [PMID 21720873](https://pubmed.ncbi.nlm.nih.gov/21720873/)
   > The HGD variant spectrum. Basis for the model's RESACT parameter (residual enzyme activity 0–5%) corresponding to the missense/null distinction.

8. Zatková A et al. (2000). *High frequency of alkaptonuria in Slovakia: evidence for the appearance of multiple mutations in HGO involving different mutational hot spots.* Am J Hum Genet. [PMID 11017803](https://pubmed.ncbi.nlm.nih.gov/11017803/)
   > The high frequency in Slovakia (about 1:19,000) and the appearance of multiple variants. Basis for the PREV/FOUNDER nodes.

9. Srsen S et al. (2002). *Alkaptonuria in Slovakia: thirty-two years of research on phenotype and genotype.* Mol Genet Metab. [PMID 12051967](https://pubmed.ncbi.nlm.nih.gov/12051967/)

10. Srsen S et al. (1978). *Screening for alkaptonuria in the newborn in Slovakia.* Lancet. [PMID 79941](https://pubmed.ncbi.nlm.nih.gov/79941/)

11. Srsen S et al. (1978). *Alkaptonuria in the Trencín District of Czechoslovakia.* Am J Med Genet. [PMID 263435](https://pubmed.ncbi.nlm.nih.gov/263435/)

12. Zatkova A et al. (2012). *Identification of 11 Novel Homogentisate 1,2 Dioxygenase Variants in Alkaptonuria Patients and Establishment of a Novel LOVD-Based HGD Mutation Database.* JIMD Rep. [PMID 23430897](https://pubmed.ncbi.nlm.nih.gov/23430897/)

13. Danda S et al. (2020). *Founder effects of the homogentisate 1,2-dioxygenase (HGD) gene in a gypsy population and mutation spectrum in the gene among alkaptonuria patients from India.* Clin Rheumatol. [PMID 32212000](https://pubmed.ncbi.nlm.nih.gov/32212000/)

14. Soltysova A et al. (2022). *Alkaptonuria in Russia.* Eur J Hum Genet. [PMID 34504318](https://pubmed.ncbi.nlm.nih.gov/34504318/)

15. Müller CR et al. (1999). *Allelic heterogeneity of alkaptonuria in Central Europe.* Eur J Hum Genet. [PMID 10482952](https://pubmed.ncbi.nlm.nih.gov/10482952/)

16. Tao L et al. (2022). *A novel mutation in the homogentisate 1,2 dioxygenase gene identified in Chinese Hani pediatric patients with Alkaptonuria.* Clin Chim Acta. [PMID 35550814](https://pubmed.ncbi.nlm.nih.gov/35550814/)

17. Abraham SSC et al. (2024). *Gene expression & biochemical analysis in alkaptonuria caused by a founder pathogenic variant across different age groups from India.* Indian J Med Res. [PMID 39737503](https://pubmed.ncbi.nlm.nih.gov/39737503/)

18. Cicaloni V et al. (2019). *Interactive alkaptonuria database: investigating clinical data to improve patient care in a rare disease.* FASEB J. [PMID 31462106](https://pubmed.ncbi.nlm.nih.gov/31462106/)


## 3. Natural history · severity scoring

19. Phornphutkul C et al. (2002). *Natural history of alkaptonuria.* N Engl J Med. [PMID 12501223](https://pubmed.ncbi.nlm.nih.gov/12501223/)
   > The source of every clinical time-point anchor in this model: joint replacement at a mean age of 55, cardiac valve involvement at 54, coronary calcium at 59, renal stones at 64, and radiographic scores rising after age 30
   > (faster in men). The model does not take these five time points in as parameters; it uses them only for validation.

20. Cox TF et al. (2011). *A quantitative assessment of alkaptonuria: testing the reliability of two disease severity scoring systems.* J Inherit Metab Dis. [PMID 21744089](https://pubmed.ncbi.nlm.nih.gov/21744089/)
   > The original AKUSSI development paper. Its three-domain clinical/joint/spine structure corresponds to the model's decomposition into the AKJ·AKS·AKC domains.

21. Langford B et al. (2018). *Alkaptonuria Severity Score Index Revisited: Analysing the AKUSSI and Its Subcomponent Features.* JIMD Rep. [PMID 29654544](https://pubmed.ncbi.nlm.nih.gov/29654544/)
   > Principal component analysis of the AKUSSI sub-items — which items carry information and which are redundant.

22. Cant HEO et al. (2022). *Improving the clinical accuracy and flexibility of the Alkaptonuria severity score index.* JIMD Rep. [PMID 35822087](https://pubmed.ncbi.nlm.nih.gov/35822087/)
   > cAKUSSI 2.0. The revision that removed biased, low-information measurements.

23. Khedr M et al. (2023). *First decade anniversary of the United Kingdom National Alkaptonuria Centre.* JIMD Rep. [PMID 36873083](https://pubmed.ncbi.nlm.nih.gov/36873083/)
   > The first 10 years of the UK National Alkaptonuria Centre (NAC). Basis for the 2 mg low-dose real-world cohort.

24. Azami A et al. (2014). *Alkaptonuric ochronosis: a clinical study from Ardabil, Iran.* Int J Rheum Dis. [PMID 24447956](https://pubmed.ncbi.nlm.nih.gov/24447956/)

25. Akbaba AI et al. (2020). *Presentation of 14 alkaptonuria patients from Turkey.* J Pediatr Endocrinol Metab. [PMID 31927521](https://pubmed.ncbi.nlm.nih.gov/31927521/)

26. Srsen S (1979). *Alkaptonuria.* Johns Hopkins Med J. [PMID 513428](https://pubmed.ncbi.nlm.nih.gov/513428/)

27. Ranganath LR et al. (2024). *Joint replacement risk is markedly increased in alkaptonuria (AKU) in those with prior arthroplasty.* Mol Genet Metab Rep. [PMID 38846518](https://pubmed.ncbi.nlm.nih.gov/38846518/)
   > The risk of a further replacement is markedly increased in patients with a previous joint replacement — consistent with the model's HZJR acceleration (cartilage loss accelerates the contralateral side through load redistribution).


## 4. Tyrosine pathway flux and mass balance — BALANCE 1 of this model

28. Milan AM et al. (2019). *Quantification of the flux of tyrosine pathway metabolites during nitisinone treatment of Alkaptonuria.* Sci Rep. [PMID 31296884](https://pubmed.ncbi.nlm.nih.gov/31296884/)
   > Quantifies the flow of tyrosine pathway metabolites using nitisinone as a flux probe. Direct basis for the model's structure in which 'the flux is conserved and only the exit changes'.

29. Ranganath LR et al. (2022). *Revisiting Quantification of Phenylalanine/Tyrosine Flux in the Ochronotic Pathway during Long-Term Nitisinone Treatment of Alkaptonuria.* Metabolites. [PMID 36295821](https://pubmed.ncbi.nlm.nih.gov/36295821/)
   > A re-examination of the paper above. Untreated, HGA equivalents amount to 7.1 g/day (46.4% of total Phe/Tyr flux); on treatment a total of 15.32 g/day is revealed, of which 8.21
   > g/day (53.6%) is unaccounted for. The authors propose biliary and intestinal excretion. This model holds that the unaccounted fraction is too large to be read as irreversible excretion (over 200 kg accumulated in 70 years), and
   > offers, quantitatively, an alternative reading in which an expanded tyrosine pool and reversible transamination are more plausible (see 'the two readings the model adjudicates between' in the README).

30. Norman BP et al. (2022). *Comprehensive Biotransformation Analysis of Phenylalanine-Tyrosine Metabolism Reveals Alternative Routes of Metabolite Clearance in Nitisinone-Treated Alkaptonuria.* Metabolites. [PMID 36295829](https://pubmed.ncbi.nlm.nih.gov/36295829/)
   > Quantification of the bypass excretion routes during treatment: urinary HPPA 20–23-fold, HPLA 9.4–10.6-fold, N-acetyl-tyrosine 5.7–6.8-fold, tyrosine glucuronide 8-fold, tyrosine
   > 4–4.7-fold. Basis for the capacity allocation across the model's tyrosine escape routes (J_TYRU·J_CONJ·J_HPPU·J_HPLAU).

31. Ranganath LR et al. (2022). *Determinants of tyrosinaemia during nitisinone therapy in alkaptonuria.* Sci Rep. [PMID 36167967](https://pubmed.ncbi.nlm.nih.gov/36167967/)
   > Analysis of 307 sampling points from SONIA 2. Concludes that what governs the degree of tyrosinaemia is not the dose but the relative decrease in the HPPA→HPLA conversion, that is, the capacity of the escape route. Independently consistent with the model's
   > asymmetry that 'the dose sets HGA, while diet and the escape routes set tyrosine'.

32. Norman BP et al. (2022). *Metabolomic studies in the inborn error of metabolism alkaptonuria reveal new biotransformations in tyrosine metabolism.* Genes Dis. [PMID 35685462](https://pubmed.ncbi.nlm.nih.gov/35685462/)

33. Norman BP et al. (2019). *A Comprehensive LC-QTOF-MS Metabolic Phenotyping Strategy: Application to Alkaptonuria.* Clin Chem. [PMID 30782595](https://pubmed.ncbi.nlm.nih.gov/30782595/)

34. Davison AS et al. (2019). *Evaluation of the serum metabolome of patients with alkaptonuria before and after two years of treatment with nitisinone using LC-QTOF-MS.* JIMD Rep. [PMID 31392115](https://pubmed.ncbi.nlm.nih.gov/31392115/)

35. Grasso D et al. (2022). *Untargeted NMR Metabolomics Reveals Alternative Biomarkers and Pathways in Alkaptonuria.* Int J Mol Sci. [PMID 36555443](https://pubmed.ncbi.nlm.nih.gov/36555443/)

36. Serafimov K et al. (2025). *Targeted and untargeted urinary metabolomics of alkaptonuria patients using ultra high-performance liquid chromatography-tandem mass spectrometry.* J Pharm Biomed Anal. [PMID 39842076](https://pubmed.ncbi.nlm.nih.gov/39842076/)

37. Gertsman I et al. (2015). *Perturbations of tyrosine metabolism promote the indolepyruvate pathway via tryptophan in host and microbiome.* Mol Genet Metab. [PMID 25680927](https://pubmed.ncbi.nlm.nih.gov/25680927/)

38. Milan AM et al. (2017). *The effect of nitisinone on homogentisic acid and tyrosine: a two-year survey of patients attending the National Alkaptonuria Centre, Liverpool.* Ann Clin Biochem. [PMID 28081634](https://pubmed.ncbi.nlm.nih.gov/28081634/)
   > Real-world HGA and tyrosine values from 2 years of NAC follow-up — the source of the 2 mg cohort anchor.


## 5. Renal handling of HGA — BALANCE 3 of this model

39. Ranganath LR et al. (2020). *Homogentisic acid is not only eliminated by glomerular filtration and tubular secretion but also produced in the kidney in alkaptonuria.* J Inherit Metab Dis. [PMID 31609457](https://pubmed.ncbi.nlm.nih.gov/31609457/)
   > An analysis of 225 patients. The renal clearance and fractional excretion of HGA exceed the theoretical maximum renal plasma flow → tubular secretion is greater than glomerular filtration and the kidney itself produces HGA. Circulating HGA also
   > rises with age, and this is significantly associated with a fall in HGA clearance. The sole basis for three terms of the model — the ceiling on CL_HGA, FREN (intrarenal production) and RFUN (age-dependent decline) —
   > and the basis for the claim that 'pigment calling for pigment is not the only reason untreated progression accelerates in old age'.

40. Wolff F et al. (2015). *Renal and prostate stones composition in alkaptonuria: a case report.* Clin Nephrol. [PMID 26396096](https://pubmed.ncbi.nlm.nih.gov/26396096/)
   > Compositional analysis of renal and prostatic stones in AKU.

41. Huledal G et al. (2019). *Non randomized study on the potential of nitisinone to inhibit cytochrome P450 2C9, 2D6, 2E1 and the organic anion transporters OAT1 and OAT3 in healthy volunteers.* Eur J Clin Pharmacol. [PMID 30443705](https://pubmed.ncbi.nlm.nih.gov/30443705/)
   > A study of the potential for nitisinone to inhibit CYP2C9/2D6/2E1 and the organic anion transporters (OAT). Basis for the pharmacological plausibility of the OAT competition term the model introduces
   > to explain the divergence between urinary and serum HGA.


## 6. Oxidation · polymerisation · pigment deposition and preclinical models

42. Preston AJ et al. (2014). *Ochronotic osteoarthropathy in a mouse model of alkaptonuria, and its inhibition by nitisinone.* Ann Rheum Dis. [PMID 23511227](https://pubmed.ncbi.nlm.nih.gov/23511227/)
   > Ochronotic osteoarthropathy in the AKU mouse model and its suppression by nitisinone. The in
   > vivo basis for the model's structure, in which only the rate of pigment deposition can be lowered while pigment already deposited cannot be reversed.

43. Hughes JH et al. (2021). *Anatomical Distribution of Ochronotic Pigment in Alkaptonuric Mice is Associated with Calcified Cartilage Chondrocytes at Osteochondral Interfaces.* Calcif Tissue Int. [PMID 33057760](https://pubmed.ncbi.nlm.nih.gov/33057760/)
   > In the mouse, the anatomical distribution of pigment is associated with the chondrocytes of calcified cartilage. Basis for the model allocating pigment preferentially to cartilage, disc and valve (FCART·FDISC·FVALV) and distinguishing these from
   > superficial tissues.

44. Ranganath LR et al. (2020). *Reversal of ochronotic pigmentation in alkaptonuria following nitisinone therapy: Analysis of data from the United Kingdom National Alkaptonuria Centre.* JIMD Rep. [PMID 32904992](https://pubmed.ncbi.nlm.nih.gov/32904992/)
   > 'Reversal' of ochronotic pigment after nitisinone treatment. This model treats it not as a refutation but as a prediction: a loss term can exist only in tissues where collagen is genuinely turned over (ear · sclera · skin)
   > and cannot exist in articular cartilage or the disc, so it is a testable consequence of the model that the observed reversal is confined to superficial tissues (in the model a
   > first-order loss term is placed only on PSCL·PEAR·PSKIN).

45. Geminiani M et al. (2017). *Cytoskeleton Aberrations in Alkaptonuric Chondrocytes.* J Cell Physiol. [PMID 27454006](https://pubmed.ncbi.nlm.nih.gov/27454006/)

46. Tinti L et al. (2011). *A novel ex vivo organotypic culture model of alkaptonuria-ochronosis.* Clin Exp Rheumatol. [PMID 21813063](https://pubmed.ncbi.nlm.nih.gov/21813063/)

47. Tinti L et al. (2010). *Evaluation of antioxidant drugs for the treatment of ochronotic alkaptonuria in an in vitro human cell model.* J Cell Physiol. [PMID 20648626](https://pubmed.ncbi.nlm.nih.gov/20648626/)

48. Braconi D et al. (2010). *Proteomic and redox-proteomic evaluation of homogentisic acid and ascorbic acid effects on human articular chondrocytes.* J Cell Biochem. [PMID 20665660](https://pubmed.ncbi.nlm.nih.gov/20665660/)
   > Proteomic and redox effects of HGA and ascorbate on human articular chondrocytes. Basis for the model explaining structurally why ascorbate reverses BQA and yet fails clinically
   > (polymerisation is a parallel, irreversible step).

49. Mastroeni P et al. (2024). *An in vitro cell model for exploring inflammatory and amyloidogenic events in alkaptonuria.* J Cell Physiol. [PMID 39351877](https://pubmed.ncbi.nlm.nih.gov/39351877/)

50. Mastroeni P et al. (2026). *HGA-Induced Oxidative Stress Impairs Autophagy via Lysosomal Dysfunction in Alkaptonuria.* Antioxidants (Basel). [PMID 42510554](https://pubmed.ncbi.nlm.nih.gov/42510554/)

51. Davison AS et al. (2022). *Impact of Nitisinone on the Cerebrospinal Fluid Metabolome of a Murine Model of Alkaptonuria.* Metabolites. [PMID 35736410](https://pubmed.ncbi.nlm.nih.gov/35736410/)

52. Davison AS et al. (2019). *Assessing the effect of nitisinone induced hypertyrosinaemia on monoamine neurotransmitters in brain tissue from a murine model of alkaptonuria using mass spectrometry imaging.* Metabolomics. [PMID 31037385](https://pubmed.ncbi.nlm.nih.gov/31037385/)
   > Nitisinone-induced hypertyrosinaemia and monoamine neurotransmitters in mouse brain tissue. Corresponds to the model's NEUROCOG node.

53. Gallagher JA et al. (2015). *Lessons from rare diseases of cartilage and bone.* Curr Opin Pharmacol. [PMID 25978274](https://pubmed.ncbi.nlm.nih.gov/25978274/)


## 7. Amyloid · inflammation · oxidative-stress biomarkers

54. Millucci L et al. (2012). *Alkaptonuria is a novel human secondary amyloidogenic disease.* Biochim Biophys Acta. [PMID 22850426](https://pubmed.ncbi.nlm.nih.gov/22850426/)
   > The proposal that AKU is a secondary amyloidogenic disease.

55. Millucci L et al. (2015). *Amyloidosis in alkaptonuria.* J Inherit Metab Dis. [PMID 25868666](https://pubmed.ncbi.nlm.nih.gov/25868666/)

56. Millucci L et al. (2014). *Amyloidosis, inflammation, and oxidative stress in the heart of an alkaptonuric patient.* Mediators Inflamm. [PMID 24876668](https://pubmed.ncbi.nlm.nih.gov/24876668/)

57. Braconi D et al. (2017). *Homogentisic acid induces aggregation and fibrillation of amyloidogenic proteins.* Biochim Biophys Acta Gen Subj. [PMID 27865997](https://pubmed.ncbi.nlm.nih.gov/27865997/)

58. Mastroeni P et al. (2024). *HGA Triggers SAA Aggregation and Accelerates Fibril Formation in the C20/A4 Alkaptonuria Cell Model.* Cells. [PMID 39273071](https://pubmed.ncbi.nlm.nih.gov/39273071/)

59. Braconi D et al. (2022). *Effects of Nitisinone on Oxidative and Inflammatory Markers in Alkaptonuria: Results from SONIA1 and SONIA2 Studies.* Cells. [PMID 36429096](https://pubmed.ncbi.nlm.nih.gov/36429096/)
   > The effect of nitisinone on oxidative and inflammatory markers in SONIA 1 and SONIA 2. Consistent with the model's structure, in which the synovitis amplification loop barely raises systemic inflammatory markers (CRP close to
   > normal).

60. Braconi D et al. (2018). *Inflammatory and oxidative stress biomarkers in alkaptonuria: data from the DevelopAKUre project.* Osteoarthritis Cartilage. [PMID 29852277](https://pubmed.ncbi.nlm.nih.gov/29852277/)

61. Trezza A et al. (2025). *Integrated Clinomics and Molecular Dynamics Simulation Approaches Reveal the SAA1.1 Allele as a Biomarker in Alkaptonuria Disease Severity.* Biomolecules. [PMID 40001497](https://pubmed.ncbi.nlm.nih.gov/40001497/)


## 8. Joint and spinal lesions and orthopaedic treatment

62. Ranganath LR et al. (2021). *Characterising the arthroplasty in spondyloarthropathy in a large cohort of eighty-seven patients with alkaptonuria.* J Inherit Metab Dis. [PMID 33314212](https://pubmed.ncbi.nlm.nih.gov/33314212/)
   > An analysis of the characteristics of joint replacement in 87 AKU patients — the principal basis for validating the timing of the model's joint replacement hazard (HZJR).

63. Salem KH et al. (2025). *Ochronotic arthropathy: skeletal manifestations and orthopaedic treatment.* EFORT Open Rev. [PMID 40071956](https://pubmed.ncbi.nlm.nih.gov/40071956/)

64. Borman P et al. (2002). *Ochronotic arthropathy.* Rheumatol Int. [PMID 11958438](https://pubmed.ncbi.nlm.nih.gov/11958438/)

65. Wu K et al. (2019). *Musculoskeletal manifestations of alkaptonuria: A case report and literature review.* Eur J Rheumatol. [PMID 30451653](https://pubmed.ncbi.nlm.nih.gov/30451653/)

66. Al-Mahfoudh R et al. (2008). *Alkaptonuria presenting with ochronotic spondyloarthropathy.* Br J Neurosurg. [PMID 19085367](https://pubmed.ncbi.nlm.nih.gov/19085367/)

67. Sang P et al. (2025). *Alkaptonuria presenting as lumbar degenerative disease: A case report and literature review.* Medicine (Baltimore). [PMID 39833049](https://pubmed.ncbi.nlm.nih.gov/39833049/)

68. Pesciallo C et al. (2022). *Total Knee Replacement in Alkaptonuric Ochronosis.* Acta Biomed. [PMID 35671127](https://pubmed.ncbi.nlm.nih.gov/35671127/)

69. Patel VG (2015). *Total knee arthroplasty in ochronosis.* Arthroplast Today. [PMID 28326376](https://pubmed.ncbi.nlm.nih.gov/28326376/)

70. Basanagoudar PL et al. (2025). *Total Joint Arthroplasty in Ochronotic Arthritis of Lower Extremities.* Indian J Orthop. [PMID 40852538](https://pubmed.ncbi.nlm.nih.gov/40852538/)

71. Tanios M et al. (2023). *Spondyloarthropathies That Mimic Ankylosing Spondylitis: A Narrative Review.* Clin Med Insights Arthritis Musculoskelet Disord. [PMID 37533960](https://pubmed.ncbi.nlm.nih.gov/37533960/)

72. Fisher AA et al. (2004). *Alkaptonuric ochronosis with aortic valve and joint replacements and femoral fracture: a case report and literature review.* Clin Med Res. [PMID 15931360](https://pubmed.ncbi.nlm.nih.gov/15931360/)


## 9. Cardiovascular involvement

73. Bruce C et al. (2025). *Effect of Nitisinone on Aortic Stenosis Disease Progression in Patients With Alkaptonuria: An Analysis of the Suitability of Nitisinone in Alkaptonuria (SONIA) 2 Study.* Cureus. [PMID 41589138](https://pubmed.ncbi.nlm.nih.gov/41589138/)
   > Analysis of aortic stenosis progression in SONIA 2. Aortic stenosis at baseline 13.0% (mild 44.4% · moderate 33.3% · severe 22.2%), aortic sclerosis 18.1%, and over 4 years the between-group difference in the rate of maximum
   > pressure-gradient progression 0.0093 mmHg/year (p=0.53, not significant). The held-out target the model's valve axis (PVALV→VALVCA→PMAXS) has to reproduce.

74. Ather N et al. (2020). *Cardiovascular ochronosis.* Cardiovasc Pathol. [PMID 32473412](https://pubmed.ncbi.nlm.nih.gov/32473412/)

75. Thimmapuram R et al. (2020). *Aortic distensibility in alkaptonuria.* Mol Genet Metab. [PMID 32466960](https://pubmed.ncbi.nlm.nih.gov/32466960/)


## 10. Ocular · aural · cutaneous signs

76. Lindner M et al. (2014). *On the ocular findings in ochronosis: a systematic review of literature.* BMC Ophthalmol. [PMID 24479547](https://pubmed.ncbi.nlm.nih.gov/24479547/)
   > A systematic review of the ocular findings in AKU — the timing and distribution of scleral pigment (Osler's sign).

77. Carlson DM et al. (1991). *Ocular ochronosis from alkaptonuria.* J Am Optom Assoc. [PMID 1813514](https://pubmed.ncbi.nlm.nih.gov/1813514/)

78. Okutucu M et al. (2019). *Glaucoma With Alkaptonuria as a Result of Pigment Accumulation.* J Glaucoma. [PMID 31274704](https://pubmed.ncbi.nlm.nih.gov/31274704/)

79. de Azevedo Magalhaes O et al. (2022). *Descemet's membrane folds in ochronosis: a case report.* J Med Case Rep. [PMID 36183119](https://pubmed.ncbi.nlm.nih.gov/36183119/)

80. Pau HW (1984). *[Involvement of the tympanic membrane and ear ossicle system in ochronotic alkaptonuria].* Laryngol Rhinol Otol (Stuttg). [PMID 6503572](https://pubmed.ncbi.nlm.nih.gov/6503572/)

81. Yin ES et al. (2017). *A 54-year-old woman with arthritis and discoloration of the hands, ears, and sclerae.* Int J Dermatol. [PMID 27651033](https://pubmed.ncbi.nlm.nih.gov/27651033/)


## 11. Bone metabolism · tendon

82. Ranganath LR et al. (2021). *Frequency, diagnosis, pathogenesis and management of osteoporosis in alkaptonuria: data analysis from the UK National Alkaptonuria Centre.* Osteoporos Int. [PMID 33118050](https://pubmed.ncbi.nlm.nih.gov/33118050/)
   > The frequency, mechanism and management of osteoporosis in AKU from the NAC data. Basis for the model's BMD term (chronic synovitis-mediated bone loss).

83. Abate M et al. (2016). *Tendons Involvement in Congenital Metabolic Disorders.* Adv Exp Med Biol. [PMID 27535253](https://pubmed.ncbi.nlm.nih.gov/27535253/)

84. Yokoe T et al. (2024). *Direct repair of the chronic ochronotic Achilles tendon rupture: a case report.* BMC Musculoskelet Disord. [PMID 39448996](https://pubmed.ncbi.nlm.nih.gov/39448996/)


## 12. HPD enzymology and the molecular pharmacology of the inhibitor

85. Ellis MK et al. (1995). *Inhibition of 4-hydroxyphenylpyruvate dioxygenase by 2-(2-nitro-4-trifluoromethylbenzoyl)-cyclohexane-1,3-dione and 2-(2-chloro-4-methanesulfonylbenzoyl)-cyclohexane-1,3-dione.* Toxicol Appl Pharmacol. [PMID 7597701](https://pubmed.ncbi.nlm.nih.gov/7597701/)
   > The original paper on the kinetics of HPD inhibition by NTBC. Basis for comparing the model's KI_NT (fitted at about 16 nmol/L) with the reported IC50 level.

86. Kavana M et al. (2003). *Interaction of (4-hydroxyphenyl)pyruvate dioxygenase with the specific inhibitor 2-[2-nitro-4-(trifluoromethyl)benzoyl]-1,3-cyclohexanedione.* Biochemistry. [PMID 12939152](https://pubmed.ncbi.nlm.nih.gov/12939152/)
   > The interaction between HPD and NTBC — slow, tight binding. Basis for the model using a dose-response steeper than simple competitive inhibition (Hill exponent > 1).

87. Lin HY et al. (2019). *Molecular insights into the mechanism of 4-hydroxyphenylpyruvate dioxygenase inhibition: enzyme kinetics, X-ray crystallography and computational simulations.* FEBS J. [PMID 30632699](https://pubmed.ncbi.nlm.nih.gov/30632699/)

88. Santucci A et al. (2017). *4-Hydroxyphenylpyruvate Dioxygenase and Its Inhibition in Plants and Animals: Small Molecules as Herbicides and Agents for the Treatment of Human Inherited Diseases.* J Med Chem. [PMID 28128559](https://pubmed.ncbi.nlm.nih.gov/28128559/)

89. Liu YX et al. (2020). *Identification of key residues determining the binding specificity of human 4-hydroxyphenylpyruvate dioxygenase.* Eur J Pharm Sci. [PMID 32750420](https://pubmed.ncbi.nlm.nih.gov/32750420/)

90. Yang DY (2003). *4-Hydroxyphenylpyruvate dioxygenase as a drug discovery target.* Drug News Perspect. [PMID 14668946](https://pubmed.ncbi.nlm.nih.gov/14668946/)

91. Brownlee JM et al. (2010). *Product analysis and inhibition studies of a causative Asn to Ser variant of 4-hydroxyphenylpyruvate dioxygenase suggest a simple route to the treatment of Hawkinsinuria.* Biochemistry. [PMID 20677779](https://pubmed.ncbi.nlm.nih.gov/20677779/)


## 13. Nitisinone pharmacokinetics and pharmacodynamics

92. Olsson B et al. (2015). *Relationship Between Serum Concentrations of Nitisinone and Its Effect on Homogentisic Acid and Tyrosine in Patients with Alkaptonuria.* JIMD Rep. [PMID 25772318](https://pubmed.ncbi.nlm.nih.gov/25772318/)
   > The relationship between serum nitisinone concentration and the effects on HGA and tyrosine. Direct basis for the model's PK (about 1 umol/L at 2 mg, about 5 umol/L at 10 mg) and for its exposure–response link.

93. Guffon N et al. (2018). *Open-Label Single-Sequence Crossover Study Evaluating Pharmacokinetics, Efficacy, and Safety of Once-Daily Dosing of Nitisinone in Patients with Hereditary Tyrosinemia Type 1.* JIMD Rep. [PMID 28643275](https://pubmed.ncbi.nlm.nih.gov/28643275/)
   > A crossover study of the PK, efficacy and safety of once-daily dosing. Basis for the model's prediction that a half-life of about 54 hours makes the dosing interval immaterial (S03 versus the divided-dose control).

94. Hall MG et al. (2001). *Pharmacokinetics and pharmacodynamics of NTBC (2-(2-nitro-4-fluoromethylbenzoyl)-1,3-cyclohexanedione) and mesotrione, inhibitors of 4-hydroxyphenyl pyruvate dioxygenase (HPPD) following a single dose to healthy male volunteers.* Br J Clin Pharmacol. [PMID 11488774](https://pubmed.ncbi.nlm.nih.gov/11488774/)

95. Anonymous (2002). *Nitisinone. Ntbc, Orfadin.* Drugs R D. [PMID 12001819](https://pubmed.ncbi.nlm.nih.gov/12001819/)

96. Gertsman I et al. (2015). *Metabolic Effects of Increasing Doses of Nitisinone in the Treatment of Alkaptonuria.* JIMD Rep. [PMID 25665838](https://pubmed.ncbi.nlm.nih.gov/25665838/)
   > The metabolic effects of raising the nitisinone dose — data showing where the dose-response saturates.

97. Kienstra NS et al. (2018). *Daily variation of NTBC and its relation to succinylacetone in tyrosinemia type 1 patients comparing a single dose to two doses a day.* J Inherit Metab Dis. [PMID 29170874](https://pubmed.ncbi.nlm.nih.gov/29170874/)

98. Schlune A et al. (2012). *Single dose NTBC-treatment of hereditary tyrosinemia type I.* J Inherit Metab Dis. [PMID 22307209](https://pubmed.ncbi.nlm.nih.gov/22307209/)


## 14. Clinical trials and regulatory status

99. Ranganath LR et al. (2016). *Suitability Of Nitisinone In Alkaptonuria 1 (SONIA 1): an international, multicentre, randomised, open-label, no-treatment controlled, parallel-group, dose-response study to investigate the effect of once daily nitisinone on 24-h urinary homogentisic acid excretion in patients with alkaptonuria after 4 weeks of treatment.* Ann Rheum Dis. [PMID 25475116](https://pubmed.ncbi.nlm.nih.gov/25475116/)
   > SONIA 1 (randomised, 5 arms, 4 weeks). Adjusted geometric mean 24-hour urinary HGA at week 4: untreated 31.53 mmol, 1 mg 3.26, 2 mg 1.44, 4 mg
   > 0.57, 8 mg 0.15 mmol. A 98.8% reduction from baseline at 8 mg. The source of all five primary anchors used to calibrate this model.

100. Ranganath LR et al. (2020). *Efficacy and safety of once-daily nitisinone for patients with alkaptonuria (SONIA 2): an international, multicentre, open-label, randomised controlled trial.* Lancet Diabetes Endocrinol. [PMID 32822600](https://pubmed.ncbi.nlm.nih.gov/32822600/)
   > SONIA 2 (randomised, 4 years, 10 mg/day, n=138). At 12 months u-HGA24 was reduced 99.7% relative to control (geometric mean ratio 0.003); at 48 months the rise in cAKUSSI
   > was 8.6 points less than in control (-16.0 ~ -1.2, p=0.023). The model's held-out clinical target.

101. Introne WJ et al. (2011). *A 3-year randomized therapeutic trial of nitisinone in alkaptonuria.* Mol Genet Metab. [PMID 21620748](https://pubmed.ncbi.nlm.nih.gov/21620748/)
   > The 3-year randomised NIH trial (2 mg/day). HGA fell substantially but the primary endpoint (hip rotation) did not improve significantly. The model has to reproduce this 'negative result',
   > and it gives the reason quantitatively (starting at age 46, most of the preventable integral is already spent).

102. Ranganath LR et al. (2022). *Comparing nitisinone 2 mg and 10 mg in the treatment of alkaptonuria-An approach using statistical modelling.* JIMD Rep. [PMID 35028273](https://pubmed.ncbi.nlm.nih.gov/35028273/)
   > A statistical modelling comparison of 2 mg and 10 mg. Serum nitisinone 1.26 vs 4.34 umol/L, urinary HGA reduced about 94% vs 99.5%, serum HGA
   > 26.2→3.86 vs 30.3→2.23 umol/L, serum tyrosine 53.9→782 vs 65.3→875 umol/L (the difference in the rate of rise between doses not significant), keratopathy
   > 3/60 (5%) vs 10/69 (14.5%), AKUSSI slope 0.19 vs 0.06 points/month. The decisive data supporting this model's central claim (the dose sets HGA, diet sets tyrosine),
   > and the object of the model's objection that the serum/urinary HGA pair at 2 mg is not compatible with mass balance.

103. Sloboda N et al. (2019). *Efficacy of low dose nitisinone in the management of alkaptonuria.* Mol Genet Metab. [PMID 31235217](https://pubmed.ncbi.nlm.nih.gov/31235217/)

104. Abbas K et al. (2022). *Adequacy of nitisinone for the management of alkaptonuria.* Ann Med Surg (Lond). [PMID 36045846](https://pubmed.ncbi.nlm.nih.gov/36045846/)

105. Ranganath LR et al. (2023). *Clinical development innovation in rare diseases: overcoming barriers to successful delivery of a randomised clinical trial in alkaptonuria-a mini-review.* Orphanet J Rare Dis. [PMID 36600285](https://pubmed.ncbi.nlm.nih.gov/36600285/)

106. Chandani HK et al. (2025). *Harliku (nitisinone): first FDA-approved disease-modifying therapy for alkaptonuria.* Ann Med Surg (Lond). [PMID 41377225](https://pubmed.ncbi.nlm.nih.gov/41377225/)
   > FDA approval of nitisinone for the AKU indication (Harliku). Following the EMA approval, the first disease-modifying treatment to be established regulatorily.

107. Ooi N et al. (2023). *Evaluation of Homogentisic Acid, a Prospective Antibacterial Agent Highlighted by the Suitability of Nitisinone in Alkaptonuria 2 (SONIA 2) Clinical Trial.* Cells. [PMID 37443717](https://pubmed.ncbi.nlm.nih.gov/37443717/)


## 15. Tyrosinaemia · keratopathy · dietary management

108. Hughes JH et al. (2020). *Dietary restriction of tyrosine and phenylalanine lowers tyrosinemia associated with nitisinone therapy of alkaptonuria.* J Inherit Metab Dis. [PMID 31503358](https://pubmed.ncbi.nlm.nih.gov/31503358/)
   > Tyrosine falls significantly with protein restriction alone, and further when Tyr/Phe-free amino acid supplementation is added (4 of 10 patients reached below 700 umol/L). Direct validation data for the model's prediction that 'diet
   > sets tyrosine', and the basis for the S15/S17 scenarios.

109. Olsson B et al. (2022). *Effects of a protein-restricted diet on body weight and serum tyrosine concentrations in patients with alkaptonuria.* JIMD Rep. [PMID 35028270](https://pubmed.ncbi.nlm.nih.gov/35028270/)
   > The effect of a protein-restricted diet on body weight and serum tyrosine.

110. Ranganath LR et al. (2024). *Anthropometric, Body Composition, and Nutritional Indicators with and without Nutritional Intervention during Nitisinone Therapy in Alkaptonuria.* Nutrients. [PMID 39203858](https://pubmed.ncbi.nlm.nih.gov/39203858/)
   > Anthropometry, body composition and nutritional indices during nitisinone treatment, with and without nutritional intervention.

111. Teke Kisa P et al. (2022). *Efficacy of Phenylalanine- and Tyrosine-Restricted Diet in Alkaptonuria Patients on Nitisinone Treatment: Case Series and Review of Literature.* Ann Nutr Metab. [PMID 34736252](https://pubmed.ncbi.nlm.nih.gov/34736252/)
   > A case series and literature review of the Phe/Tyr-restricted diet. The source of the actual protocol — 0.9 g/kg at a tyrosine of 501–700, 0.8 g/kg at 701–900, and the addition of Phe/Tyr-free exchanges
   > above 900 — and the basis for setting the dietary levels of the model's S14–S17.

112. Ranganath LR et al. (2022). *Comparing the Phenylalanine/Tyrosine Pathway and Related Factors between Keratopathy and No-Keratopathy Groups as Well as between Genders in Alkaptonuria during Nitisinone Treatment.* Metabolites. [PMID 36005644](https://pubmed.ncbi.nlm.nih.gov/36005644/)
   > A comparison of the Phe/Tyr pathway between the keratopathy and non-keratopathy groups. The data at which the threshold structure of the model's CORTYR (corneal crystal burden) is aimed.

113. Ranganath L et al. (2023). *Increased prevalence of Parkinson's disease in alkaptonuria.* JIMD Rep. [PMID 37404676](https://pubmed.ncbi.nlm.nih.gov/37404676/)
   > An increased prevalence of Parkinson's disease in AKU. A potential link through the tyrosine-dopamine axis (the model's NEUROCOG node, marked uncertain).

114. Rudebeck M et al. (2020). *A patient survey on the impact of alkaptonuria symptoms as perceived by the patients and their experiences of receiving diagnosis and care.* JIMD Rep. [PMID 32395411](https://pubmed.ncbi.nlm.nih.gov/32395411/)

115. Davison AS et al. (2016). *Acute fatal metabolic complications in alkaptonuria.* J Inherit Metab Dis. [PMID 26596578](https://pubmed.ncbi.nlm.nih.gov/26596578/)
   > An acute fatal metabolic complication of AKU — cited in order to state explicitly that this is a route the model does not cover, rare though it is.


## 16. Long-term nitisinone experience in hereditary tyrosinaemia type 1 (basis for extrapolation)

116. Spiekerkoetter U et al. (2021). *Long-term safety and outcomes in hereditary tyrosinaemia type 1 with nitisinone treatment: a 15-year non-interventional, multicentre study.* Lancet Diabetes Endocrinol. [PMID 34023005](https://pubmed.ncbi.nlm.nih.gov/34023005/)
   > A 15-year non-interventional long-term safety study. The safety profile observed under lifelong exposure from childhood — the only human evidence the model's S09 (start at age 5) scenario leans on, and the reason it must be stated that this extrapolation
   > has not been tested in AKU.

117. Geppert J et al. (2017). *Evaluation of pre-symptomatic nitisinone treatment on long-term outcomes in Tyrosinemia type 1 patients: a systematic review.* Orphanet J Rare Dis. [PMID 28893311](https://pubmed.ncbi.nlm.nih.gov/28893311/)
   > Long-term outcomes of nitisinone treatment started before symptom onset — the evidence from HT-1 on the value of early initiation.

118. Chinsky JM et al. (2017). *Diagnosis and treatment of tyrosinemia type I: a US and Canadian consensus group review and recommendations.* Genet Med. [PMID 28771246](https://pubmed.ncbi.nlm.nih.gov/28771246/)

119. van Ginkel WG et al. (2019). *Long-Term Outcomes and Practical Considerations in the Pharmacological Management of Tyrosinemia Type 1.* Paediatr Drugs. [PMID 31667718](https://pubmed.ncbi.nlm.nih.gov/31667718/)

120. Kuypers AM et al. (2025). *Overview of European Practices for Management of Tyrosinemia Type 1: Towards European Guidelines.* J Inherit Metab Dis. [PMID 40965374](https://pubmed.ncbi.nlm.nih.gov/40965374/)

121. Holme E et al. (1995). *Diagnosis and management of tyrosinemia type I.* Curr Opin Pediatr. [PMID 8776026](https://pubmed.ncbi.nlm.nih.gov/8776026/)

122. Aktuglu Zeybek AC et al. (2022). *Evaluation of dynamic thiol/disulfide homeostasis in hereditary tyrosinemia type 1 patients.* Pediatr Res. [PMID 34628487](https://pubmed.ncbi.nlm.nih.gov/34628487/)


## 17. Ascorbate and other attempted therapies

123. Wolff JA et al. (1989). *Effects of ascorbic acid in alkaptonuria: alterations in benzoquinone acetic acid and an ontogenic effect in infancy.* Pediatr Res. [PMID 2771520](https://pubmed.ncbi.nlm.nih.gov/2771520/)
   > Ascorbate alters BQA but does not change the clinical course; an age effect in infancy. Basis for the model reproducing the failure of ascorbate structurally (irreversible polymerisation is a parallel route).

124. Mayatepek E et al. (1998). *Effects of ascorbic acid and low-protein diet in alkaptonuria.* Eur J Pediatr. [PMID 9809834](https://pubmed.ncbi.nlm.nih.gov/9809834/)

125. Morava E et al. (2003). *Reversal of clinical symptoms and radiographic abnormalities with protein restriction and ascorbic acid in alkaptonuria.* Ann Clin Biochem. [PMID 12542920](https://pubmed.ncbi.nlm.nih.gov/12542920/)

126. Forslind K et al. (1988). *Alkaptonuria and ochronosis in three siblings. Ascorbic acid treatment monitored by urinary HGA excretion.* Clin Exp Rheumatol. [PMID 3180550](https://pubmed.ncbi.nlm.nih.gov/3180550/)

127. Kamoun P et al. (1992). *Ascorbic acid and alkaptonuria.* Eur J Pediatr. [PMID 1537362](https://pubmed.ncbi.nlm.nih.gov/1537362/)


## 18. Investigational treatment strategies

128. Lequeue S et al. (2025). *A robust bacterial high-throughput screening assay to identify pharmacological chaperones targeting human homogentisate 1,2-dioxygenase missense variants in alkaptonuria.* Eur J Pharmacol. [PMID 40784658](https://pubmed.ncbi.nlm.nih.gov/40784658/)
   > A bacterium-based high-throughput screen to find pharmacological chaperones targeting HGD. Corresponds to the model's CHAPERONE→RESACT route, and shares the structural constraint that it is meaningful only in
   > missense genotypes.

129. Paulk NK et al. (2012). *In vivo selection of transplanted hepatocytes by pharmacological inhibition of fumarylacetoacetate hydrolase in wild-type mice.* Mol Ther. [PMID 22871666](https://pubmed.ncbi.nlm.nih.gov/22871666/)

130. Vernon HJ et al. (2021). *Milestones in treatments for inborn errors of metabolism: Reflections on Where chemistry and medicine meet.* Am J Med Genet A. [PMID 34165242](https://pubmed.ncbi.nlm.nih.gov/34165242/)


## 19. QSP methodology

131. Elmokadem A et al. (2019). *Quantitative Systems Pharmacology and Physiologically-Based Pharmacokinetic Modeling With mrgsolve: A Hands-On Tutorial.* CPT Pharmacometrics Syst Pharmacol. [PMID 31652028](https://pubmed.ncbi.nlm.nih.gov/31652028/)
   > A hands-on tutorial on building QSP and PBPK models with mrgsolve — the tooling basis for this model's implementation.

---

## Where no citation is claimed

The values below are **not** cited directly from the literature; they were derived inside the model or set
arbitrarily. They are collected here so that they do not appear to have an evidential basis.

- **Pigment partition fractions** (`FCART` 0.30 · `FDISC` 0.22 · `FVALV` 0.08 · `FTEND` 0.12 ·
  `FSCL` 0.03 · `FEAR` 0.05 · `FSKIN` 0.05 · `FOTH` 0.15): allocated from the qualitative frequencies of
  involvement in the autopsy and imaging descriptive literature, preserving their rank order only; there is no quantitative basis.
- **Collagen binding capacities** (`COLLC` 3000 · `COLLD` 2200 · `COLLV` 800 · `COLLT` 1500 umol):
  scaling constants that normalise pigment density onto a 0–1 range. The absolute values carry no meaning; they
  act only as a single ratio together with `PD50`.
- **Embrittlement Hill exponent** `HPD_H` = 6: chosen as the least steep slope able to reproduce the
  observation that symptoms appear abruptly in the third decade. There is no histological basis.
- **Embrittlement thresholds** `PD50` 0.60 / `PD50D` 0.45 / `PD50V` 0.70 / `PD50T` 0.80:
  an ordering constraint based only on the clinical sequence in which disc calcification appears first and
  tendon rupture last; the individual values are arbitrary.
- **AKUSSI domain weights** (45 / 40 / 30 / 15): a reduction of the real AKUSSI's 57 sub-item structure to
  4 domains, which does not reproduce the item-by-item scoring of the original scale. The model's absolute
  cAKUSSI values therefore cannot be compared directly with clinical scores; **only the rate of change (points/month)**
  is comparable. This is why the validation table uses slopes rather than absolute values.
- **Pain module coefficients** (`WNOCI_*`, `KCS`, `KCSOFF`): with no pain psychophysics data specific to
  AKU, only the general form from the chronic osteoarthritis literature was borrowed.
- **Nitisinone peripheral compartment** (`V2`, `Q`): values that keep a two-compartment structure, not
  estimated from a population PK analysis in AKU patients. The steady-state concentration is set by `CLNT` alone.
- **`SRC_INS` = 7.5 umol/day**: a source of HGA that does not respond to nitisinone was assumed and fitted,
  and the fit **converged on almost zero**. That is, the data support the hypothesis that 'the dose–response
  is steeper than simple competitive inhibition' (`HNT` 1.45) over the 'uninhibitable HGA source' hypothesis.
  It is left in as a term that was hypothesised and then rejected by the data.

---

## Documented disagreements between the literature and the model

1. **Serum HGA at 2 mg.** Reported value 3.86 umol/L. The model gives 15.1 umol/L (3.9-fold too high).
   This is not a parameter problem but a mass-balance problem. Essentially all HGA leaves in the urine,
   so urinary HGA ≈ total HGA production. If urinary output is 1,200–1,440 umol/day at 2 mg and
   181 umol/day at 10 mg, then for serum to be 3.86 and 2.23 the HGA clearance would have to be
   **4.6-fold higher** at 2 mg than at 10 mg. Yet the concentrations of the competing organic anions (HPPA·HPLA)
   are almost the same at the two doses. No single clearance model can satisfy the four values
   simultaneously. Candidate explanations are (i) the limit of quantification of serum HGA at low concentrations,
   (ii) a cohort difference between NAC (2 mg, active dietary management) and SONIA 2 (10 mg, information only),
   and (iii) a saturable secretory route not yet described. The discriminating experiment: **measure serum and urinary
   HGA at 2 mg and at 10 mg simultaneously, in the same cohort at the same visit** (addressed to PMID 35028273).
2. **Magnitude of the serum HPPA rise.** Reported 14.3–15.0-fold, model 24.9-fold. The steep inhibition
   curve at `HNT` = 1.45 pushes HPP up further than it should. The same 1.70-fold overestimate appears in the
   2 mg/10 mg ratio of urinary HGA (reported 6.63, model 11.3), and the two errors being the same multiple is
   not a coincidence but a consequence of both being set by the same steepness of inhibition.
3. **The effect of nitisinone on the aortic valve.** SONIA 2 reported no effect, with a between-group difference in
   the rate of maximum pressure-gradient progression over 4 years of 0.0093 mmHg/year (p=0.53). The model excludes the valve
   route from the reversible channel and gives -0.057 mmHg/year — the sign is right and the magnitude still 6-fold
   too high, but set against the control group's 4-year progression (about 4 mmHg) it is 1.4%, effectively null.
4. **The rate of reversal of pigment deposition.** PMID 32904992 reports a visible reversal of skin, ear and
   scleral pigment after nitisinone treatment. The model reduces it by only 0.5% over 2 years (expected about 20%). This is because the model's
   dermal turnover half-life of 6 years cannot produce a visible change in only 2 years. The observed
   reversal implies that **the effective turnover of visible pigment is on a scale of months**, and this is a measurable
   parameter (serial photographic measurement + the rate of re-deposition after nitisinone is stopped).
5. **The timing of renal stone formation.** Urinary HGA supersaturation is present from birth, so under a first-order
   reaction model stones form in infancy. Gating on renal papillary pigment as the nucleus improved this to
   age 59.4 (reported 64), but the model has no structure for a stone nucleation induction time,
   so 'the point at which it becomes clinically apparent' depends on an arbitrary detection threshold.

---

## Tools

- **mrgsolve** — <https://mrgsolve.org> · <https://github.com/metrumresearchgroup/mrgsolve>
- **Graphviz** — <https://graphviz.org>
- **R / Shiny** — <https://shiny.posit.co>
- **PubMed E-utilities** — <https://www.ncbi.nlm.nih.gov/books/NBK25501/>
