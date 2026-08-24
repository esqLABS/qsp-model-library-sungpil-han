# Central Precocious Puberty QSP Model — References
# Central Precocious Puberty (CPP) QSP Model — References

This document collects **the grounds for every structural assumption and parameter** that went into `cpp_qsp_model.dot`, `cpp_mrgsolve_model_en.R`, `cpp_shiny_app_en.R` and
`cpp_reference_check.py`.

**97 papers in total** · every PMID was looked up through the PubMed E-utilities and its authors, journal, year and title cross-checked
(the bibliographic information returned by the verification script was used exactly as returned, and no PMID here was written from memory).

The note on each entry records "which **parameter, or which structural choice**, of the model
this paper supports". Literature that is merely cited and not used in the model has not been included.

---

## 1. Genetics of pubertal onset — CPP begins with a LOST BRAKE

Human puberty is held in a state of active inhibition, and most monogenic CPP is not a gain of drive but a **loss of inhibition**. In the model MKRN3(0) = `MK0` is the only patient-level knob, and the age of pubertal onset is an output.

1. **Abreu AP** et al. *Central precocious puberty caused by mutations in the imprinted gene MKRN3.* N Engl J Med. 2013.  
   [PMID 23738509](https://pubmed.ncbi.nlm.nih.gov/23738509/)  
   → The first report that loss of MKRN3 function is the commonest monogenic cause of familial CPP. The reason the model's `MK0` parameter exists.

2. **Macedo DB** et al. *Central precocious puberty that appears to be sporadic caused by paternally inherited mutations in the imprinted gene makorin ring finger 3.* J Clin Endocrinol Metab. 2014.  
   [PMID 24628548](https://pubmed.ncbi.nlm.nih.gov/24628548/)  
   → A paternally inherited MKRN3 variant is found even in CPP that looks 'sporadic' — because the gene is imprinted, maternal transmission produces no phenotype.

3. **Macedo DB** et al. *Central Precocious Puberty Caused by a Heterozygous Deletion in the MKRN3 Promoter Region.* Neuroendocrinology. 2018.  
   [PMID 29763903](https://pubmed.ncbi.nlm.nih.gov/29763903/)  
   → A deletion in the MKRN3 promoter region — the brake can be lost without a coding variant.

4. **Stecchini MF** et al. *Time Course of Central Precocious Puberty Development Caused by an MKRN3 Gene Mutation: A Prismatic Case.* Horm Res Paediatr. 2016.  
   [PMID 27424312](https://pubmed.ncbi.nlm.nih.gov/27424312/)  
   → A case following the time course of MKRN3-mutant CPP — the grounds for `TAUKND` of 60 days.

5. **Abreu AP** et al. *A new pathway in the control of the initiation of puberty: the MKRN3 gene.* J Mol Endocrinol. 2015.  
   [PMID 25957321](https://pubmed.ncbi.nlm.nih.gov/25957321/)  
   → A review of the MKRN3 pathway.

6. **Busch AS** et al. *Circulating MKRN3 Levels Decline During Puberty in Healthy Boys.* J Clin Endocrinol Metab. 2016.  
   [PMID 27057785](https://pubmed.ncbi.nlm.nih.gov/27057785/)  
   → Shows that circulating MKRN3 falls before puberty in healthy boys — the model's `KMK` exponential decay (half-life 3.0 years) comes from here.

7. **Roberts SA** et al. *Hypothalamic Overexpression of Makorin Ring Finger Protein 3 Results in Delayed Puberty in Female Mice.* Endocrinology. 2022.  
   [PMID 35974456](https://pubmed.ncbi.nlm.nih.gov/35974456/)  
   → Hypothalamic MKRN3 overexpression delays puberty in the female mouse — the experimental grounds for the model's direction of effect, `MK0`↑ → delayed puberty.

8. **Dauber A** et al. *Paternally Inherited DLK1 Deletion Associated With Familial Central Precocious Puberty.* J Clin Endocrinol Metab. 2017.  
   [PMID 28324015](https://pubmed.ncbi.nlm.nih.gov/28324015/)  
   → A paternal DLK1 deletion causes familial CPP plus central obesity — the model's `DLK1R` residual inhibition term.

9. **Teles MG** et al. *A GPR54-activating mutation in a patient with central precocious puberty.* N Engl J Med. 2008.  
   [PMID 18272894](https://pubmed.ncbi.nlm.nih.gov/18272894/)  
   → CPP from an activating KISS1R(GPR54) variant — the counterexample, a gain-of-drive form.

10. **Silveira LG** et al. *Mutations of the KISS1 gene in disorders of puberty.* J Clin Endocrinol Metab. 2010.  
   [PMID 20237166](https://pubmed.ncbi.nlm.nih.gov/20237166/)  
   → KISS1 gene variants and disorders of puberty.

11. **de Roux N** et al. *Hypogonadotropic hypogonadism due to loss of function of the KiSS1-derived peptide receptor GPR54.* Proc Natl Acad Sci U S A. 2003.  
   [PMID 12944565](https://pubmed.ncbi.nlm.nih.gov/12944565/)  
   → Loss of GPR54 function → hypogonadotropic hypogonadism. The opposite extreme of the same axis.

12. **Topaloglu AK** et al. *TAC3 and TACR3 mutations in familial hypogonadotropic hypogonadism reveal a key role for Neurokinin B in the central control of reproduction.* Nat Genet. 2009.  
   [PMID 19079066](https://pubmed.ncbi.nlm.nih.gov/19079066/)  
   → TAC3/TACR3 variants — establish that NKB is essential to KNDy self-excitation (the model's pulse ignition).

13. **Canton APM** et al. *Rare variants in the MECP2 gene in girls with central precocious puberty: a translational cohort study.* Lancet Diabetes Endocrinol. 2023.  
   [PMID 37385287](https://pubmed.ncbi.nlm.nih.gov/37385287/)  
   → Rare MECP2 variants and CPP in girls.

14. **Canton APM** et al. *MECP2 Rare Variants in Boys With Central Precocious Puberty.* J Clin Endocrinol Metab. 2026.  
   [PMID 41139199](https://pubmed.ncbi.nlm.nih.gov/41139199/)  
   → Rare MECP2 variants and CPP in boys.

15. **Day FR** et al. *Genomic analyses identify hundreds of variants associated with age at menarche and support a role for puberty timing in cancer risk.* Nat Genet. 2017.  
   [PMID 28436984](https://pubmed.ncbi.nlm.nih.gov/28436984/)  
   → A GWAS of hundreds of variants associated with age at menarche — the genetic share of the variance in pubertal timing.

16. **Abreu AP** et al. *Pubertal development and regulation.* Lancet Diabetes Endocrinol. 2016.  
   [PMID 26852256](https://pubmed.ncbi.nlm.nih.gov/26852256/)  
   → Pubertal development and its regulation — a review covering the whole axis.

17. **Roberts SA** et al. *GENETICS IN ENDOCRINOLOGY: Genetic etiologies of central precocious puberty and the role of imprinted genes.* Eur J Endocrinol. 2020.  
   [PMID 32698138](https://pubmed.ncbi.nlm.nih.gov/32698138/)  
   → A review of the genetic causes of CPP and of imprinted genes.

---

## 2. The KNDy pulse generator — the state variable is not concentration but FREQUENCY

The model holds the **frequency** of GnRH pulses (`PULS`) as an explicit state variable. NKB ignites a pulse and dynorphin terminates it, and the interval between pulses follows from that.

18. **Goodman RL** et al. *Kisspeptin neurons in the arcuate nucleus of the ewe express both dynorphin A and neurokinin B.* Endocrinology. 2007.  
   [PMID 17823266](https://pubmed.ncbi.nlm.nih.gov/17823266/)  
   → Arcuate kisspeptin neurons in the sheep co-express dynorphin A and NKB — the anatomical grounds for the KNDy concept.

19. **Merkley CM** et al. *KNDy (kisspeptin/neurokinin B/dynorphin) neurons are activated during both pulsatile and surge secretion of LH in the ewe.* Endocrinology. 2012.  
   [PMID 22989631](https://pubmed.ncbi.nlm.nih.gov/22989631/)  
   → KNDy neurons are activated both in pulsatile secretion and in the LH surge.

20. **Clarkson J** et al. *Definition of the hypothalamic GnRH pulse generator in mice.* Proc Natl Acad Sci U S A. 2017.  
   [PMID 29109258](https://pubmed.ncbi.nlm.nih.gov/29109258/)  
   → Optogenetic definition of the hypothalamic GnRH pulse generator in the mouse — proof that KNDy population activity is itself the pulse.

21. **McQuillan HJ** et al. *Definition of the estrogen negative feedback pathway controlling the GnRH pulse generator in female mice.* Nat Commun. 2022.  
   [PMID 36460649](https://pubmed.ncbi.nlm.nih.gov/36460649/)  
   → Definition of the oestrogen negative feedback route that controls the pulse generator — the model's `KFB` term.

22. **Garcia JP** et al. *Kisspeptin and Neurokinin B Signaling Network Underlies the Pubertal Increase in GnRH Release in Female Rhesus Monkeys.* Endocrinology. 2017.  
   [PMID 28977601](https://pubmed.ncbi.nlm.nih.gov/28977601/)  
   → Shows that in the primate the kisspeptin-NKB signalling network underlies the pubertal rise in GnRH secretion.

23. **Toro CA** et al. *Trithorax dependent changes in chromatin landscape at enhancer and promoter regions drive female puberty.* Nat Commun. 2018.  
   [PMID 29302059](https://pubmed.ncbi.nlm.nih.gov/29302059/)  
   → Trithorax-dependent chromatin change drives female puberty — the epigenetic switch.

24. **Vazquez MJ** et al. *SIRT1 mediates obesity- and nutrient-dependent perturbation of pubertal timing by epigenetically controlling Kiss1 expression.* Nat Commun. 2018.  
   [PMID 30305620](https://pubmed.ncbi.nlm.nih.gov/30305620/)  
   → SIRT1 mediates the nutrition/obesity-dependent disturbance of pubertal timing — the model's `KLEPB` leptin permissive signal.

25. **Bessa DS** et al. *Methylome profiling of healthy and central precocious puberty girls.* Clin Epigenetics. 2018.  
   [PMID 30466473](https://pubmed.ncbi.nlm.nih.gov/30466473/)  
   → Methylome profiling of girls with CPP.

26. **Teilmann G** et al. *Increased risk of precocious puberty in internationally adopted children in Denmark.* Pediatrics. 2006.  
   [PMID 16882780](https://pubmed.ncbi.nlm.nih.gov/16882780/)  
   → The increased risk of precocious puberty in internationally adopted children — the epidemiology behind the model's metabolic/environmental permissive signal.

---

## 3. Pulsatility decoding and GnRHR pharmacology — an agonist does not block, it 'destroys the code'

The literature in this section is what justifies the model's central equation `S = RS·(Sendo·ffree + AINT·fa)`. Because of its intrinsic activity `AINT` = 1.6, an agonist at first **raises** stimulation (the flare), and suppression comes only after the sensitised receptor pool `RS` has collapsed.

27. **Belchetz PE** et al. *Hypophysial responses to continuous and intermittent delivery of hypopthalamic gonadotropin-releasing hormone.* Science. 1978.  
   [PMID 100883](https://pubmed.ncbi.nlm.nih.gov/100883/)  
   → The classic experiment from Knobil's group: continuous administration suppresses the gonadotropins while pulsatile administration sustains them. The therapeutic principle of the entire model is this one paper.

28. **Wildt L** et al. *Frequency and amplitude of gonadotropin-releasing hormone stimulation and gonadotropin secretion in the rhesus monkey.* Endocrinology. 1981.  
   [PMID 6788538](https://pubmed.ncbi.nlm.nih.gov/6788538/)  
   → In the rhesus monkey the frequency and amplitude of GnRH stimulation determine gonadotropin secretion — the source of the `PULSMIN`/`PULSMAX` range.

29. **Conn PM** et al. *Gonadotropin-releasing hormone and its analogues.* N Engl J Med. 1991.  
   [PMID 1984190](https://pubmed.ncbi.nlm.nih.gov/1984190/)  
   → GnRH and its analogues — a review of the clinical pharmacology.

30. **Millar RP** et al. *Gonadotropin-releasing hormone receptors.* Endocr Rev. 2004.  
   [PMID 15082521](https://pubmed.ncbi.nlm.nih.gov/15082521/)  
   → A review of the GnRH receptor. That the type II receptor carries no C-terminal tail, so that recycling is slow and incomplete, is the basis for `KREC` = 0.035/day.

31. **Sealfon SC** et al. *Molecular mechanisms of ligand interaction with the gonadotropin-releasing hormone receptor.* Endocr Rev. 1997.  
   [PMID 9101136](https://pubmed.ncbi.nlm.nih.gov/9101136/)  
   → The molecular mechanism of the GnRH receptor-ligand interaction — the grounds for the `EC50` values.

32. **Kaiser UB** et al. *Differential effects of gonadotropin-releasing hormone (GnRH) pulse frequency on gonadotropin subunit and GnRH receptor messenger ribonucleic acid levels in vitro.* Endocrinology. 1997.  
   [PMID 9048630](https://pubmed.ncbi.nlm.nih.gov/9048630/)  
   → The differential effect of GnRH pulse frequency on gonadotropin subunit and receptor mRNA — fast pulses give LHB, slow pulses FSHB. The model's pulsatility decoder.

33. **Kanasaki H** et al. *Gonadotropin-releasing hormone pulse frequency-dependent activation of extracellular signal-regulated kinase pathways in perifused LbetaT2 cells.* Endocrinology. 2005.  
   [PMID 16141398](https://pubmed.ncbi.nlm.nih.gov/16141398/)  
   → GnRH pulse frequency-dependent activation of the ERK pathway.

34. **Tornøe CW** et al. *Pharmacokinetic/pharmacodynamic modelling of GnRH antagonist degarelix: a comparison of the non-linear mixed-effects programs NONMEM and NLME.* J Pharmacokinet Pharmacodyn. 2004.  
   [PMID 16222784](https://pubmed.ncbi.nlm.nih.gov/16222784/)  
   → PK/PD modelling of the GnRH antagonist degarelix — the structural reference for the model's antagonist arm (`KAX`, `EC50X`).

35. **Tornøe CW** et al. *Population pharmacokinetic/pharmacodynamic (PK/PD) modelling of the hypothalamic-pituitary-gonadal axis following treatment with GnRH analogues.* Br J Clin Pharmacol. 2007.  
   [PMID 17096678](https://pubmed.ncbi.nlm.nih.gov/17096678/)  
   → Population PK/PD modelling of the hypothalamic-pituitary-gonadal axis after degarelix — the quantitative reference for the kinetics of suppression.

---

## 4. The growth plate — only oestrogen closes it, and it does so not by time but by a DEPLETED reserve

The (−) arm of the model. Human loss-of-function experiments (ESR1, aromatase deficiency) established that the end of skeletal maturation is **oestrogen-dependent** in both sexes, and Baron's group formalised that as the depletion of a reserve rather than as a clock. The model's `GPRES` is exactly this.

36. **Morishima A** et al. *Aromatase deficiency in male and female siblings caused by a novel mutation and the physiological role of estrogens.* J Clin Endocrinol Metab. 1995.  
   [PMID 8530621](https://pubmed.ncbi.nlm.nih.gov/8530621/)  
   → Aromatase deficiency in a brother and a sister — without oestrogen the growth plate does not close.

37. **Carani C** et al. *Effect of testosterone and estradiol in a man with aromatase deficiency.* N Engl J Med. 1997.  
   [PMID 9211678](https://pubmed.ncbi.nlm.nih.gov/9211678/)  
   → A comparison of the effects of testosterone and of oestradiol in an aromatase-deficient man — what closes skeletal maturation is E2, not T. The grounds on which the model may use an aromatase inhibitor in boys as a 'tool for separating the two signs'.

38. **Nilsson O** et al. *Fundamental limits on longitudinal bone growth: growth plate senescence and epiphyseal fusion.* Trends Endocrinol Metab. 2004.  
   [PMID 15380808](https://pubmed.ncbi.nlm.nih.gov/15380808/)  
   → The fundamental limit on longitudinal growth: growth plate senescence and the end of skeletal maturation. The conceptual source of `GPRES` 1 → 0.

39. **Weise M** et al. *Effects of estrogen on growth plate senescence and epiphyseal fusion.* Proc Natl Acad Sci U S A. 2001.  
   [PMID 11381135](https://pubmed.ncbi.nlm.nih.gov/11381135/)  
   → Experimental demonstration that oestrogen accelerates growth plate senescence and the end of skeletal maturation.

40. **Nilsson O** et al. *Evidence that estrogen hastens epiphyseal fusion and cessation of longitudinal bone growth by irreversibly depleting the number of resting zone progenitor cells in female rabbits.* Endocrinology. 2014.  
   [PMID 24708243](https://pubmed.ncbi.nlm.nih.gov/24708243/)  
   → Direct evidence that oestrogen brings skeletal maturation forward by **irreversibly depleting** the growth plate progenitor reserve — why `GPRES` never recovers in the model.

41. **Emons J** et al. *Mechanisms of growth plate maturation and epiphyseal fusion.* Horm Res Paediatr. 2011.  
   [PMID 21540578](https://pubmed.ncbi.nlm.nih.gov/21540578/)  
   → A review of the mechanisms of growth plate maturation and of the end of skeletal maturation.

42. **Grumbach MM** et al. *Estrogen, bone, growth and sex: a sea change in conventional wisdom.* J Pediatr Endocrinol Metab. 2000.  
   [PMID 11202221](https://pubmed.ncbi.nlm.nih.gov/11202221/)  
   → Oestrogen · bone · growth — a commentary setting out the turn in received opinion.

43. **Emons J** et al. *Expression of vascular endothelial growth factor in the growth plate is stimulated by estradiol and increases during pubertal development.* J Endocrinol. 2010.  
   [PMID 20093283](https://pubmed.ncbi.nlm.nih.gov/20093283/)  
   → Growth plate VEGF expression is stimulated by oestradiol and rises at puberty.

44. **BAYLEY N** et al. *Tables for predicting adult height from skeletal age: revised for use with the Greulich-Pyle hand standards.* J Pediatr. 1952.  
   [PMID 14918032](https://pubmed.ncbi.nlm.nih.gov/14918032/)  
   → The Bayley-Pinneau tables — the predicted adult height that is calculated in the clinic. The model keeps these tables separate, so that the 'clinical prediction' and 'the final height the model actually integrated' sit side by side.

---

## 5. Guidelines and final adult-height outcomes

The reference literature against which the model's sign-flip result is checked. The consensus statement records that there is no height gain when treatment starts after a bone age of 12, and the model reproduces that point **by computing** it, at a bone age of 11.6 years (a gain of under 1 cm).

45. **Carel JC** et al. *Consensus statement on the use of gonadotropin-releasing hormone analogs in children.* Pediatrics. 2009.  
   [PMID 19332438](https://pubmed.ncbi.nlm.nih.gov/19332438/)  
   → The international consensus statement on the use of GnRH analogues in children — the reference for treatment indications, monitoring, and when to stop.

46. **Carel JC** et al. *Clinical practice. Precocious puberty.* N Engl J Med. 2008.  
   [PMID 18509122](https://pubmed.ncbi.nlm.nih.gov/18509122/)  
   → A clinical review of precocious puberty.

47. **Latronico AC** et al. *Causes, diagnosis, and treatment of central precocious puberty.* Lancet Diabetes Endocrinol. 2016.  
   [PMID 26852255](https://pubmed.ncbi.nlm.nih.gov/26852255/)  
   → The causes, diagnosis and treatment of central precocious puberty.

48. **Bangalore Krishna K** et al. *Use of Gonadotropin-Releasing Hormone Analogs in Children: Update by an International Consortium.* Horm Res Paediatr. 2019.  
   [PMID 31319416](https://pubmed.ncbi.nlm.nih.gov/31319416/)  
   → The international consortium update on the use of GnRH analogues.

49. **Latronico AC** et al. *Central precocious puberty: an Endocrine Society clinical practice guideline.* J Clin Endocrinol Metab. 2026.  
   [PMID 42287186](https://pubmed.ncbi.nlm.nih.gov/42287186/)  
   → Central precocious puberty — the Endocrine Society clinical practice guideline (the most recent).

50. **Carel JC** et al. *Final height after long-term treatment with triptorelin slow release for central precocious puberty: importance of statural growth after interruption of treatment. French study group of Decapeptyl in Precocious Puberty.* J Clin Endocrinol Metab. 1999.  
   [PMID 10372696](https://pubmed.ncbi.nlm.nih.gov/10372696/)  
   → Final height after long-term treatment with sustained-release triptorelin — the principal reference value for final height in the treated group.

51. **Klein KO** et al. *Increased final height in precocious puberty after long-term treatment with LHRH agonists: the National Institutes of Health experience.* J Clin Endocrinol Metab. 2001.  
   [PMID 11600530](https://pubmed.ncbi.nlm.nih.gov/11600530/)  
   → The gain in final height after long-term LHRH agonist treatment: the NIH experience.

52. **Heger S** et al. *Long-term outcome after depot gonadotropin-releasing hormone agonist treatment of central precocious puberty: final height, body proportions, body composition, bone mineral density, and reproductive function.* J Clin Endocrinol Metab. 1999.  
   [PMID 10599723](https://pubmed.ncbi.nlm.nih.gov/10599723/)  
   → Long-term outcomes after depot GnRH agonist treatment — final height, body proportions, body composition, bone mineral density.

53. **Lazar L** et al. *Growth pattern and final height after cessation of gonadotropin-suppressive therapy in girls with central sexual precocity.* J Clin Endocrinol Metab. 2007.  
   [PMID 17579199](https://pubmed.ncbi.nlm.nih.gov/17579199/)  
   → Growth pattern and final height after gonadal suppression is stopped — the check on the model's post-treatment recovery and timing of menarche.

54. **Pasquino AM** et al. *Long-term observation of 87 girls with idiopathic central precocious puberty treated with gonadotropin-releasing hormone analogs: impact on adult height, body mass index, bone mineral content, and reproductive function.* J Clin Endocrinol Metab. 2008.  
   [PMID 17940112](https://pubmed.ncbi.nlm.nih.gov/17940112/)  
   → Long-term observation of 87 girls with idiopathic CPP — final height, BMI, ovarian function.

55. **Antoniazzi F** et al. *End results in central precocious puberty with GnRH analog treatment: the data of the Italian Study Group for Physiopathology of Puberty.* J Pediatr Endocrinol Metab. 2000.  
   [PMID 10969920](https://pubmed.ncbi.nlm.nih.gov/10969920/)  
   → The Italian study group's final results for GnRH analogue treatment of CPP.

56. **Arrigo T** et al. *When to stop GnRH analog therapy: the experience of the Italian Study Group for Physiopathology of Puberty.* J Pediatr Endocrinol Metab. 2000.  
   [PMID 10969918](https://pubmed.ncbi.nlm.nih.gov/10969918/)  
   → When to stop a GnRH analogue — the Italian study group.

57. **Bereket A** et al. *A Critical Appraisal of the Effect of Gonadotropin-Releasing Hormon Analog Treatment on Adult Height of Girls with Central Precocious Puberty.* J Clin Res Pediatr Endocrinol. 2017.  
   [PMID 29280737](https://pubmed.ncbi.nlm.nih.gov/29280737/)  
   → A critical appraisal of the effect of GnRH analogue treatment on final height — it emphasises that the size of the gain depends strongly on when treatment starts.

58. **Guaraldi F** et al. *MANAGEMENT OF ENDOCRINE DISEASE: Long-term outcomes of the treatment of central precocious puberty.* Eur J Endocrinol. 2016.  
   [PMID 26466612](https://pubmed.ncbi.nlm.nih.gov/26466612/)  
   → A review of the long-term outcomes of treating central precocious puberty.

59. **Magiakou MA** et al. *The efficacy and safety of gonadotropin-releasing hormone analog treatment in childhood and adolescence: a single center, long-term follow-up study.* J Clin Endocrinol Metab. 2010.  
   [PMID 19897682](https://pubmed.ncbi.nlm.nih.gov/19897682/)  
   → The efficacy and safety of GnRH analogue treatment in childhood and adolescence — long-term follow-up at a single centre.

60. **Hwangbo J** et al. *Long-term outcomes of gonadotropin-releasing hormone agonist treatment in girls with central precocious puberty.* Ann Pediatr Endocrinol Metab. 2025.  
   [PMID 40049673](https://pubmed.ncbi.nlm.nih.gov/40049673/)  
   → Long-term outcomes of GnRH agonist treatment in girls with CPP (a recent cohort).

61. **Cho KW** et al. *Final adult height in male patients with central precocious puberty after gonadotropin-releasing hormone agonist treatment.* Ann Pediatr Endocrinol Metab. 2026.  
   [PMID 41787710](https://pubmed.ncbi.nlm.nih.gov/41787710/)  
   → Final adult height after GnRH agonist treatment in boys with CPP.

---

## 6. Formulation · dose · PK — the design variable is not potency but trough coverage

The model computes that there is almost no difference between depot formulations given exactly on time, and that failure arises from (a) a nasal spray, (b) a delayed injection, (c) a low dose.

62. **Eugster EA** et al. *Efficacy and safety of histrelin subdermal implant in children with central precocious puberty: a multicenter trial.* J Clin Endocrinol Metab. 2007.  
   [PMID 17327379](https://pubmed.ncbi.nlm.nih.gov/17327379/)  
   → The multicentre trial of the efficacy and safety of the histrelin subcutaneous implant — zero-order release of 65 µg/day, a formulation with no trough.

63. **Silverman LA** et al. *Long-Term Continuous Suppression With Once-Yearly Histrelin Subcutaneous Implants for the Treatment of Central Precocious Puberty: A Final Report of a Phase 3 Multicenter Trial.* J Clin Endocrinol Metab. 2015.  
   [PMID 25803268](https://pubmed.ncbi.nlm.nih.gov/25803268/)  
   → Sustained long-term suppression with the once-yearly histrelin implant.

64. **Fuld K** et al. *A randomized trial of 1- and 3-month depot leuprolide doses in the treatment of central precocious puberty.* J Pediatr. 2011.  
   [PMID 21798557](https://pubmed.ncbi.nlm.nih.gov/21798557/)  
   → A randomised comparison of 1-month against 3-month depot leuprolide — the test of the model's prediction that it effectively cannot tell the two formulations apart.

65. **Klein K** et al. *Efficacy and safety of triptorelin 6-month formulation in patients with central precocious puberty.* J Pediatr Endocrinol Metab. 2016.  
   [PMID 26887034](https://pubmed.ncbi.nlm.nih.gov/26887034/)  
   → The efficacy and safety of the 6-month triptorelin formulation.

66. **Mostafa NM** et al. *Pharmacokinetic and exposure-response analyses of leuprolide following administration of leuprolide acetate 3-month depot formulations to children with central precocious puberty.* Clin Drug Investig. 2014.  
   [PMID 24756362](https://pubmed.ncbi.nlm.nih.gov/24756362/)  
   → PK and exposure-response analysis of the 3-month leuprolide depot — the source for calibrating `FBURSTL`, `KDISL` and `CLL`.

67. **Lim CN** et al. *Fixed Dosing of Leuprolide Acetate, a GnRH Agonist, in Children with Central Precocious Puberty: A Population Pharmacokinetic Justification.* Paediatr Drugs. 2026.  
   [PMID 41824265](https://pubmed.ncbi.nlm.nih.gov/41824265/)  
   → A population PK analysis of fixed-dose leuprolide in paediatric CPP.

68. **Wang X** et al. *Reflections on GnRH analogues dosage in the context of escape phenomenon: a case report of hypothalamic hamartoma with central precocious puberty and literature review.* J Pediatr Endocrinol Metab. 2026.  
   [PMID 42381421](https://pubmed.ncbi.nlm.nih.gov/42381421/)  
   → Reconsidering the GnRH analogue dose in the context of the escape phenomenon — a case of hypothalamic hamartoma. The model's `KNDLES` ectopic drive term.

---

## 7. Diagnosis and on-treatment monitoring

What the model's monitoring tab sets out to check. **Growth velocity actually falls under effective suppression, so it cannot be trusted**, while ΔBA/ΔCA, basal LH and uterine volume do work.

69. **Neely EK** et al. *Normal ranges for immunochemiluminometric gonadotropin assays.* J Pediatr. 1995.  
   [PMID 7608809](https://pubmed.ncbi.nlm.nih.gov/7608809/)  
   → Normal ranges for the gonadotropins by immunochemiluminescence — the anchor for the model's LH calibration (prepubertal 0.1-0.15 IU/L).

70. **Houk CP** et al. *Adequacy of a single unstimulated luteinizing hormone level to diagnose central precocious puberty in girls.* Pediatrics. 2009.  
   [PMID 19482738](https://pubmed.ncbi.nlm.nih.gov/19482738/)  
   → Can CPP in a girl be diagnosed from a single unstimulated LH — specificity is high, sensitivity low.

71. **Resende EA** et al. *Assessment of basal and gonadotropin-releasing hormone-stimulated gonadotropins by immunochemiluminometric and immunofluorometric assays in normal children.* J Clin Endocrinol Metab. 2007.  
   [PMID 17284632](https://pubmed.ncbi.nlm.nih.gov/17284632/)  
   → A comparison of assays for basal and GnRH-stimulated gonadotropins — the grounds for the post-stimulation peak LH > 5 IU/L criterion.

72. **de Vries L** et al. *Role of pelvic ultrasound in girls with precocious puberty.* Horm Res Paediatr. 2011.  
   [PMID 21228561](https://pubmed.ncbi.nlm.nih.gov/21228561/)  
   → The role of pelvic ultrasound in girls with precocious puberty — the uterine volume/length criteria (> 3.4 mL / > 34 mm).

73. **Chalumeau M** et al. *Selecting girls with precocious puberty for brain imaging: validation of European evidence-based diagnosis rule.* J Pediatr. 2003.  
   [PMID 14571217](https://pubmed.ncbi.nlm.nih.gov/14571217/)  
   → Selecting which girls with precocious puberty to image the brain in — validation of the European evidence-based diagnostic rule.

74. **Cho MH** et al. *Diagnostic Usefulness of Subcutaneous Triptorelin Stimulation Testing in the Evaluation of Central Precocious Puberty in Girls.* Horm Metab Res. 2026.  
   [PMID 41956123](https://pubmed.ncbi.nlm.nih.gov/41956123/)  
   → The diagnostic usefulness of the subcutaneous triptorelin stimulation test.

---

## 8. Bone mineral density · body composition · long-term safety

The model computes that the 'fall' in bone mineral density Z-score is mostly **an advanced baseline regressing to the same-age reference**, and that peak bone mass recovers.

75. **Boot AM** et al. *Bone mineral density and body composition before and during treatment with gonadotropin-releasing hormone agonist in children with central precocious and early puberty.* J Clin Endocrinol Metab. 1998.  
   [PMID 9467543](https://pubmed.ncbi.nlm.nih.gov/9467543/)  
   → Bone mineral density and body composition before and during GnRH agonist treatment in children with CPP — the source for calibrating `KZ` and `KCATCH`.

76. **van der Sluis IM** et al. *Longitudinal follow-up of bone density and body composition in children with precocious or early puberty before, during and after cessation of GnRH agonist therapy.* J Clin Endocrinol Metab. 2002.  
   [PMID 11836277](https://pubmed.ncbi.nlm.nih.gov/11836277/)  
   → Longitudinal follow-up of bone mineral density and body composition in children with precocious or early puberty.

77. **Antoniazzi F** et al. *Prevention of bone demineralization by calcium supplementation in precocious puberty during gonadotropin-releasing hormone agonist treatment.* J Clin Endocrinol Metab. 1999.  
   [PMID 10372699](https://pubmed.ncbi.nlm.nih.gov/10372699/)  
   → Calcium supplementation prevents bone demineralisation during GnRH agonist treatment — the model's `CAVD`/`KCAVD` term.

78. **Heger S** et al. *Long-term GnRH agonist treatment for female central precocious puberty does not impair reproductive function.* Mol Cell Endocrinol. 2006.  
   [PMID 16757104](https://pubmed.ncbi.nlm.nih.gov/16757104/)  
   → Long-term GnRH agonist treatment of CPP in girls does not impair reproductive function.

79. **Franceschi R** et al. *Prevalence of polycystic ovary syndrome in young women who had idiopathic central precocious puberty.* Fertil Steril. 2010.  
   [PMID 19135667](https://pubmed.ncbi.nlm.nih.gov/19135667/)  
   → The prevalence of polycystic ovary syndrome in young women who had idiopathic CPP.

80. **Schoelwer MJ** et al. *Psychological assessment of mothers and their daughters at the time of diagnosis of precocious puberty.* Int J Pediatr Endocrinol. 2015.  
   [PMID 25780366](https://pubmed.ncbi.nlm.nih.gov/25780366/)  
   → Psychological assessment of mothers and daughters at the time of a diagnosis of precocious puberty — the grounds for the model's psychosocial index (`QOL`) term.

---

## 9. GnRH-INDEPENDENT (peripheral) precocious puberty — the territory an agonist structurally cannot reach

One of the model's sharpest results: the same leuprolide regimen gives +9.8 cm in central disease and **+0.1 cm** in McCune-Albright. Where the target does not match the mechanism, the drug effect is zero.

81. **Weinstein LS** et al. *Activating mutations of the stimulatory G protein in the McCune-Albright syndrome.* N Engl J Med. 1991.  
   [PMID 1944469](https://pubmed.ncbi.nlm.nih.gov/1944469/)  
   → The activating Gsα mutation of McCune-Albright syndrome — the mechanism by which cAMP is switched on without LH. The model's `AUTSET` term.

82. **Shenker A** et al. *Severe endocrine and nonendocrine manifestations of the McCune-Albright syndrome associated with activating mutations of stimulatory G protein GS.* J Pediatr. 1993.  
   [PMID 8410501](https://pubmed.ncbi.nlm.nih.gov/8410501/)  
   → The severe endocrine and non-endocrine phenotypes of McCune-Albright associated with an activating Gsα mutation.

83. **Eugster EA** et al. *Tamoxifen treatment for precocious puberty in McCune-Albright syndrome: a multicenter trial.* J Pediatr. 2003.  
   [PMID 12915825](https://pubmed.ncbi.nlm.nih.gov/12915825/)  
   → The multicentre trial of tamoxifen for precocious puberty in McCune-Albright — the model's `TAMEFF` term.

84. **Feuillan P** et al. *Letrozole treatment of precocious puberty in girls with the McCune-Albright syndrome: a pilot study.* J Clin Endocrinol Metab. 2007.  
   [PMID 17405850](https://pubmed.ncbi.nlm.nih.gov/17405850/)  
   → Letrozole treatment of precocious puberty in girls with McCune-Albright — the grounds for a high-potency AI giving +8.8 cm in the model.

85. **Sims EK** et al. *Fulvestrant treatment of precocious puberty in girls with McCune-Albright syndrome.* Int J Pediatr Endocrinol. 2012.  
   [PMID 22999294](https://pubmed.ncbi.nlm.nih.gov/22999294/)  
   → Fulvestrant treatment in girls with McCune-Albright.

86. **Reiter EO** et al. *Bicalutamide plus anastrozole for the treatment of gonadotropin-independent precocious puberty in boys with testotoxicosis: a phase II, open-label pilot study (BATT).* J Pediatr Endocrinol Metab. 2010.  
   [PMID 21158211](https://pubmed.ncbi.nlm.nih.gov/21158211/)  
   → Bicalutamide plus anastrozole for gonadotropin-independent precocious puberty in boys with testotoxicosis.

---

## 10. Aromatase inhibitors and rhGH — the attempt to separate the two signs

An AI blocks E2 while leaving testosterone, so in principle it removes the (−) arm alone. Yet the model computes that an AI on its own does not improve height in boys, and this agrees with the actual trial results — because taking E2 away also takes away the GH/IGF-1 amplification, that is, the (+) arm.

87. **Hero M** et al. *Inhibition of estrogen biosynthesis with a potent aromatase inhibitor increases predicted adult height in boys with idiopathic short stature: a randomized controlled trial.* J Clin Endocrinol Metab. 2005.  
   [PMID 16189252](https://pubmed.ncbi.nlm.nih.gov/16189252/)  
   → Potent aromatase inhibition raises predicted adult height in boys with idiopathic short stature.

88. **Wickman S** et al. *A specific aromatase inhibitor and potential increase in adult height in boys with delayed puberty: a randomised controlled trial.* Lancet. 2001.  
   [PMID 11403810](https://pubmed.ncbi.nlm.nih.gov/11403810/)  
   → An aromatase inhibitor and adult height in boys with delayed puberty — a randomised trial.

89. **Varimo T** et al. *Letrozole Monotherapy in Pre- and Early-Pubertal Boys Does Not Increase Adult Height.* Front Endocrinol (Lausanne). 2019.  
   [PMID 31024444](https://pubmed.ncbi.nlm.nih.gov/31024444/)  
   → **Letrozole monotherapy does not increase adult height in prepubertal and early pubertal boys** — the improvement in predicted height did not carry through into actual height. A prediction confirmed after the fact, agreeing in direction with the model's computed −0.3 cm for an AI alone.

90. **Mauras N** et al. *Anastrozole increases predicted adult height of short adolescent males treated with growth hormone: a randomized, placebo-controlled, multicenter trial for one to three years.* J Clin Endocrinol Metab. 2008.  
   [PMID 18165285](https://pubmed.ncbi.nlm.nih.gov/18165285/)  
   → Anastrozole raises predicted adult height in short boys on growth hormone treatment.

91. **Pasquino AM** et al. *Adult height in girls with central precocious puberty treated with gonadotropin-releasing hormone analogues and growth hormone.* J Clin Endocrinol Metab. 1999.  
   [PMID 10022399](https://pubmed.ncbi.nlm.nih.gov/10022399/)  
   → Adult height with a GnRH analogue plus growth hormone in girls with CPP — the reference for the size of the added rhGH gain (a few cm).

92. **Pasquino AM** et al. *Adult height in short normal girls treated with gonadotropin-releasing hormone analogs and growth hormone.* J Clin Endocrinol Metab. 2000.  
   [PMID 10690865](https://pubmed.ncbi.nlm.nih.gov/10690865/)  
   → A GnRH analogue plus growth hormone in short girls who are otherwise normal.

---

## 11. Epidemiology and secular trends

The literature that settles which population the model's normal reference child (thelarche at 10.3 years, menarche at 13.0) stands for.

93. **Teilmann G** et al. *Prevalence and incidence of precocious pubertal development in Denmark: an epidemiologic study based on national registries.* Pediatrics. 2005.  
   [PMID 16322154](https://pubmed.ncbi.nlm.nih.gov/16322154/)  
   → The prevalence and incidence of precocious puberty in Denmark — from the national registers.

94. **Soriano-Guillén L** et al. *Central precocious puberty in children living in Spain: incidence, prevalence, and influence of adoption and immigration.* J Clin Endocrinol Metab. 2010.  
   [PMID 20554707](https://pubmed.ncbi.nlm.nih.gov/20554707/)  
   → The incidence and prevalence of CPP in children resident in Spain, and the influence of adoption and immigration.

95. **Bräuner EV** et al. *Trends in the Incidence of Central Precocious Puberty and Normal Variant Puberty Among Children in Denmark, 1998 to 2017.* JAMA Netw Open. 2020.  
   [PMID 33044548](https://pubmed.ncbi.nlm.nih.gov/33044548/)  
   → Incidence trends for central precocious puberty and for normal variant puberty in Denmark, 1998-2017.

96. **Eckert-Lind C** et al. *Worldwide Secular Trends in Age at Pubertal Onset Assessed by Breast Development Among Girls: A Systematic Review and Meta-analysis.* JAMA Pediatr. 2020.  
   [PMID 32040143](https://pubmed.ncbi.nlm.nih.gov/32040143/)  
   → Worldwide secular trends in the age at onset of thelarche in girls — a systematic review and meta-analysis. The source of the model's `KMK` calibration target (thelarche at 10.3 years).

97. **Stagi S** et al. *Increased incidence of precocious and accelerated puberty in females during and after the Italian lockdown for the coronavirus 2019 (COVID-19) pandemic.* Ital J Pediatr. 2020.  
   [PMID 33148304](https://pubmed.ncbi.nlm.nih.gov/33148304/)  
   → The rise in precocious puberty and accelerated puberty in girls during and after the Italian COVID-19 lockdown.

---

## Explicitly NOT verified against the literature

Stated for the sake of honesty. The items below were not derived directly from the literature above.

- **`ENDOBLEED` = 0.58** — the **one and only** parameter in the whole model that was fitted to a published figure.
  It was set so that the incidence of withdrawal bleeding after the first depot falls inside the reported 5-10% range.
- **The psychosocial index (`QOL`) and the hot flush index (`HF`)** — composite indices on a 0-10 scale, corresponding to no
  validated instrument whatever. Only their direction and their relative magnitude should be read.
- **`KMKSEX` = 0.75 (the brake comes loose later in boys)** — a phenomenological scale for reproducing the observation
  that puberty in boys is about 1.8 years later than in girls; it is not based on any paper reporting a sex difference
  in MKRN3 itself.
- **The GnRH antagonist regimen (18 mg q28d)** — no antagonist formulation is approved for CPP. The model's antagonist
  arm is a **hypothetical formulation** for computing the counterfactual of "suppression without a flare", and
  its PK borrows the structure of the degarelix literature (PMID 16222784, 17096678).
- **`XTRA` = 0.10** — the assumption that the reserve is depleted slightly faster than the maturity a bone-age atlas
  reads. The direction agrees with PMID 24708243 but there is no established value for the magnitude.

---

## Where this model DISAGREES with the literature

Set out quantitatively in the last section of `README_en.md`. In summary there are three.

1. Peak growth velocity in the normal child is about 1 cm/year lower than in the population data (the growth plate
   reserve index, once it is made to match final height, flattens the growth spurt).
2. With early-start treatment ΔBA/ΔCA is about 0.95-1.00, higher than the commonly quoted 0.4-0.7.
   Where treatment started late it falls to 0.56-0.72, inside the published range.
3. The Bayley-Pinneau predicted adult height comes out 5-10 cm above the final height the model actually integrated.
   This is not a bug but a prediction of the model, and it runs in the same direction as the clinical observation
   that the BP method overestimates in a child whose bone age is advanced.

