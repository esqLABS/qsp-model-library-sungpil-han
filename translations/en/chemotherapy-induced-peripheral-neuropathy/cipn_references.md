# Chemotherapy-Induced Peripheral Neuropathy (CIPN) — References
# Chemotherapy-Induced Peripheral Neuropathy — Annotated References

Each entry notes which structure or parameter of the QSP model it supports.
Each entry notes which model structure or parameter it supports.

---

## 1. Epidemiology and clinical overview

1. Seretny M, et al. **Incidence, prevalence, and predictors of chemotherapy-induced peripheral neuropathy: a systematic review and meta-analysis.** *Pain.* 2014;155(12):2461-2470.
   <https://pubmed.ncbi.nlm.nih.gov/25261162/>
   — a prevalence of 68.1% (1 month) → 60% (3 months) → 30% (beyond 6 months). The basis for calibrating the recovery time constants of the model (§8 recovery).

2. Staff NP, Grisold A, Grisold W, Windebank AJ. **Chemotherapy-induced peripheral neuropathy: a current review.** *Ann Neurol.* 2017;81(6):772-781.
   <https://pubmed.ncbi.nlm.nih.gov/28486769/>
   — an overview of the mechanisms by class. The skeleton of the structure of map clusters 4, 5, and 6.

3. Colvin LA. **Chemotherapy-induced peripheral neuropathy: where are we now?** *Pain.* 2019;160(Suppl 1):S1-S10.
   <https://pubmed.ncbi.nlm.nih.gov/31008843/>
   — the distinction between the acute and chronic phenotypes, and the role of central sensitisation.

4. Loprinzi CL, et al. **Prevention and Management of Chemotherapy-Induced Peripheral Neuropathy in Survivors of Adult Cancers: ASCO Guideline Update.** *J Clin Oncol.* 2020;38(28):3325-3348.
   <https://pubmed.ncbi.nlm.nih.gov/32663120/>
   — duloxetine is the only recommended treatment. The basis for the list of agents that failed in prevention (cluster 13).

5. Bonhof CS, et al. **Course of chemotherapy-induced peripheral neuropathy and its impact on health-related quality of life among colorectal cancer patients.** *Int J Cancer.* 2019;145(9):2482-2491.
   <https://pubmed.ncbi.nlm.nih.gov/30721528/>
   — the CIPN20 trajectory and its link to QoL (cluster 15).

6. Selvy M, et al. **Long-Term Prevalence of Sensory Chemotherapy-Induced Peripheral Neuropathy for 5 Years after Adjuvant FOLFOX Chemotherapy.** *J Clin Med.* 2020;9(8):2400.
   <https://pubmed.ncbi.nlm.nih.gov/32727063/>
   — residual CIPN at 5 years. The basis for the floor of irreversible neuronal loss (the NEURON compartment).

7. Kerckhove N, et al. **Long-Term Effects, Pathophysiological Mechanisms, and Risk Factors of Chemotherapy-Induced Peripheral Neuropathies: A Comprehensive Literature Review.** *Front Pharmacol.* 2017;8:86.
   <https://pubmed.ncbi.nlm.nih.gov/28286483/>

---

## 2. DRG vulnerability and transporters

8. Jimenez-Andrade JM, et al. **Vascularization of the dorsal root ganglia and peripheral nerve of the mouse: implications for chemical-induced peripheral sensory neuropathies.** *Mol Pain.* 2008;4:10.
   <https://pubmed.ncbi.nlm.nih.gov/18353190/>
   — the permeability of the DRG vasculature (the absence of a blood-nerve barrier). The key basis for cluster 3.

9. Sprowl JA, et al. **Oxaliplatin-induced neurotoxicity is dependent on the organic cation transporter OCT2.** *Proc Natl Acad Sci USA.* 2013;110(27):11199-11204.
   <https://pubmed.ncbi.nlm.nih.gov/23776246/>
   — OCT2-mediated uptake of oxaliplatin into the DRG. The model's `OCT2` covariate.

10. Ip V, et al. **Quantification of platinum-DNA adducts in dorsal root ganglia and peripheral blood in a rat model.** *Anticancer Res.* 2012;32(3):995-1002.
    <https://pubmed.ncbi.nlm.nih.gov/22399623/>
    — quantification of platinum-DNA adducts in the DRG. The `ADDUCT` compartment.

11. Ceresa C, Cavaletti G. **Drug transporters in chemotherapy induced peripheral neurotoxicity: current knowledge and clinical implications.** *Curr Med Chem.* 2011;18(3):329-341.
    <https://pubmed.ncbi.nlm.nih.gov/21174632/>

12. Cavaletti G, et al. **Distribution of paclitaxel within the nervous system of the rat after repeated intravenous administration.** *Neurotoxicology.* 2000;21(3):389-393.
    <https://pubmed.ncbi.nlm.nih.gov/10894128/>
    — the DRG-predominant distribution of the taxanes. The `CE_TAX` effect compartment.

---

## 3. Platinum: adducts and nucleolar stress

13. McDonald ES, Randon KR, Knight A, Windebank AJ. **Cisplatin preferentially binds to DNA in dorsal root ganglion neurons in vitro and in vivo: a X-ray fluorescence microscopy study.** *Neurobiol Dis.* 2005;18(2):305-313.
    <https://pubmed.ncbi.nlm.nih.gov/15686960/>

14. Gill JS, Windebank AJ. **Cisplatin-induced apoptosis in rat dorsal root ganglion neurons is associated with attempted entry into the cell cycle.** *J Clin Invest.* 1998;101(12):2842-2850.
    <https://pubmed.ncbi.nlm.nih.gov/9637718/>
    — p53-mediated death in a non-dividing neuron. The `NEURON` loss equation.

15. Podratz JL, et al. **Cisplatin induced mitochondrial DNA damage in dorsal root ganglion neurons.** *Neurobiol Dis.* 2011;41(3):661-668.
    <https://pubmed.ncbi.nlm.nih.gov/21145397/>
    — mtDNA has no NER, so damage accumulates. The mechanistic basis for damage to the mitochondrial compartment.

16. Melli G, Taiana M, Camozzi F, et al. **Alpha-lipoic acid prevents mitochondrial damage and neurotoxicity in experimental chemotherapy neuropathy.** *Exp Neurol.* 2008;214(2):276-284.
    <https://pubmed.ncbi.nlm.nih.gov/18809400/>

17. Bruna J, et al. **Neurophysiological, histological and immunohistochemical characterization of bortezomib-induced neuropathy in mice.** *Exp Neurol.* 2010;223(2):599-608.
    <https://pubmed.ncbi.nlm.nih.gov/20188093/>

18. Screnci D, et al. **Relationships between hydrophobicity, reactivity, accumulation and peripheral nerve toxicity of a series of platinum drugs.** *Br J Cancer.* 2000;82(4):966-972.
    <https://pubmed.ncbi.nlm.nih.gov/10732772/>
    — the amount of platinum accumulated within the nerve determines the neurotoxicity. The `CE_PT` accumulation structure.

---

## 4. Microtubules and axonal transport

19. Gornstein E, Schwarz TL. **The paradox of paclitaxel neurotoxicity: Mechanisms and unanswered questions.** *Neuropharmacology.* 2014;76 Pt A:175-183.
    <https://pubmed.ncbi.nlm.nih.gov/23978385/>
    — the transport hypothesis of taxane neurotoxicity. The key basis for the `ATRANS` compartment.

20. Bobylev I, et al. **Paclitaxel inhibits mRNA transport in axons.** *Neurobiol Dis.* 2015;82:321-331.
    <https://pubmed.ncbi.nlm.nih.gov/26183709/>

21. LaPointe NE, et al. **Effects of eribulin, vincristine, paclitaxel and ixabepilone on fast axonal transport and kinesin-1 driven microtubule gliding: implications for chemotherapy-induced peripheral neuropathy.** *Neurotoxicology.* 2013;37:231-239.
    <https://pubmed.ncbi.nlm.nih.gov/23711742/>
    — differences in the strength of transport inhibition between drugs. The basis for the `KDAM` scaler by class.

22. Smith JA, et al. **Tubulin-dependent changes in the L-type calcium channel: implications for CIPN.** *Neurotoxicology.* 2016;. (review)
    <https://pubmed.ncbi.nlm.nih.gov/27036088/>

23. Pease-Raissi SE, et al. **Paclitaxel Reduces Axonal Bclw to Initiate IP3R1-Dependent Axon Degeneration.** *Neuron.* 2017;96(2):373-386.e6.
    <https://pubmed.ncbi.nlm.nih.gov/29024662/>

---

## 5. The SARM1 / NAD⁺ axon-death programme

24. Gerdts J, Brace EJ, Sasaki Y, DiAntonio A, Milbrandt J. **SARM1 activation triggers axon degeneration locally via NAD+ destruction.** *Science.* 2015;348(6233):453-457.
    <https://pubmed.ncbi.nlm.nih.gov/25908823/>
    — the mechanistic basis for the `SARM` compartment.

25. Geisler S, et al. **Prevention of vincristine-induced peripheral neuropathy by genetic deletion of SARM1 in mice.** *Brain.* 2016;139(Pt 12):3092-3108.
    <https://pubmed.ncbi.nlm.nih.gov/27797810/>
    — SARM1 deficiency prevents CIPN. The executing programme of progression after dosing (coasting).

26. Geisler S, et al. **Vincristine and bortezomib use distinct upstream mechanisms to activate a common SARM1-dependent axon degeneration program.** *JCI Insight.* 2019;4(17):e129920.
    <https://pubmed.ncbi.nlm.nih.gov/31484833/>
    — different classes converge on a common SARM1 pathway — the basis for the model using a single `SARM` compartment.

27. Gilley J, Coleman MP. **Endogenous Nmnat2 is an essential survival factor for maintenance of healthy axons.** *PLoS Biol.* 2010;8(1):e1000300.
    <https://pubmed.ncbi.nlm.nih.gov/20126265/>
    — the dependence of NMNAT2 on axonal transport → the `ATRANS`→`ENERGY`→`SARM` link.

28. Bosanac T, et al. **Pharmacological SARM1 inhibition protects axon structure and function in paclitaxel-induced peripheral neuropathy.** *Brain.* 2021;144(10):3226-3238.
    <https://pubmed.ncbi.nlm.nih.gov/33964142/>
    — an administrable SARM1 inhibitor. The investigational node of cluster 13.

29. Coleman MP, Höke A. **Programmed axon degeneration: from mouse to mechanism to medicine.** *Nat Rev Neurosci.* 2020;21(4):183-196.
    <https://pubmed.ncbi.nlm.nih.gov/32152523/>

---

## 6. The coasting phenomenon (worsening after treatment ends)

30. Argyriou AA, et al. **Either called 'chemobrain' or 'chemofog', the long-term chemotherapy-induced cognitive decline in cancer survivors is real.** — a related review that describes coasting: Argyriou AA, et al. **Chemotherapy-induced peripheral neurotoxicity (CIPN): an update.** *Crit Rev Oncol Hematol.* 2012;82(1):51-77.
    <https://pubmed.ncbi.nlm.nih.gov/21908200/>

31. Briani C, et al. **Long-term course of oxaliplatin-induced polyneuropathy: a prospective 2-year follow-up study.** *J Peripher Nerv Syst.* 2014;19(4):299-306.
    <https://pubmed.ncbi.nlm.nih.gov/25583063/>
    — partial recovery after the post-treatment worsening. The basis for calibrating `KS_OFF` (t½ 23 d) and `KREGEN_GATE`.

32. Pachman DR, et al. **Clinical Course of Oxaliplatin-Induced Neuropathy: Results From the Randomized Phase III Trial N08CB (Alliance).** *J Clin Oncol.* 2015;33(30):3416-3422.
    <https://pubmed.ncbi.nlm.nih.gov/26282635/>
    — prospective quantification of coasting and the separation of the acute from the chronic symptoms.

33. Bennedsgaard K, et al. **Oxaliplatin- and docetaxel-induced polyneuropathy: clinical and neurophysiological characteristics.** *J Peripher Nerv Syst.* 2020;25(4):377-387.
    <https://pubmed.ncbi.nlm.nih.gov/32902895/>

---

## 7. Ion channels and the acute oxaliplatin syndrome

34. Sittl R, et al. **Anticancer drug oxaliplatin induces acute cooling-aggravated neuropathy via sodium channel subtype Na(V)1.6-resurgent and persistent current.** *Proc Natl Acad Sci USA.* 2012;109(17):6704-6709.
    <https://pubmed.ncbi.nlm.nih.gov/22493249/>
    — the Nav1.6 persistent current. The molecular basis for the `COLDA` state.

35. Descoeur J, et al. **Oxaliplatin-induced cold hypersensitivity is due to remodelling of ion channel expression in nociceptors.** *EMBO Mol Med.* 2011;3(5):266-278.
    <https://pubmed.ncbi.nlm.nih.gov/21438154/>
    — TRPM8↑, Kv reduced. Cluster 9.

36. Grolleau F, et al. **A possible explanation for a neurotoxic effect of the anticancer agent oxaliplatin on neuronal voltage-gated sodium channels.** *J Neurophysiol.* 2001;85(5):2293-2297.
    <https://pubmed.ncbi.nlm.nih.gov/11353042/>
    — the Ca²⁺ chelation hypothesis of oxalate.

37. Park SB, et al. **Oxaliplatin-induced neurotoxicity: changes in axonal excitability precede development of neuropathy.** *Brain.* 2009;132(Pt 10):2712-2723.
    <https://pubmed.ncbi.nlm.nih.gov/19745023/>
    — the change in excitability precedes the neuropathy. The model structure in which `EXCITC` comes before the structural damage.

38. Deuis JR, et al. **An animal model of oxaliplatin-induced cold allodynia reveals a crucial role for Nav1.6 in peripheral pain pathways.** *Pain.* 2013;154(9):1749-1757.
    <https://pubmed.ncbi.nlm.nih.gov/23711482/>

---

## 8. Neuroinflammation

39. Zhang H, Li Y, de Carvalho-Barbosa M, et al. **Dorsal Root Ganglion Infiltration by Macrophages Contributes to Paclitaxel Chemotherapy-Induced Peripheral Neuropathy.** *J Pain.* 2016;17(7):775-786.
    <https://pubmed.ncbi.nlm.nih.gov/27063781/>
    — the `MAC` compartment.

40. Makker PGS, et al. **Characterisation of Immune and Neuroinflammatory Changes Associated with Chemotherapy-Induced Peripheral Neuropathy.** *PLoS One.* 2017;12(1):e0170814.
    <https://pubmed.ncbi.nlm.nih.gov/28125674/>

41. Li Y, et al. **Toll-Like Receptor 4 Signaling Contributes to Paclitaxel-Induced Peripheral Neuropathy.** *J Pain.* 2014;15(7):712-725.
    <https://pubmed.ncbi.nlm.nih.gov/24755282/>
    — TLR4/MyD88 → IL-1β. Cluster 10.

42. Luo X, et al. **Macrophage Toll-like Receptor 9 Contributes to Chemotherapy-Induced Neuropathic Pain in Male Mice.** *J Neurosci.* 2019;39(35):6848-6864.
    <https://pubmed.ncbi.nlm.nih.gov/31270160/>

43. Old EA, et al. **Monocytes expressing CX3CR1 orchestrate the development of vincristine-induced pain.** *J Clin Invest.* 2014;124(5):2023-2036.
    <https://pubmed.ncbi.nlm.nih.gov/24743151/>

---

## 9. Central sensitisation and descending modulation

44. Boyette-Davis JA, Walters ET, Dougherty PM. **Mechanisms involved in the development of chemotherapy-induced neuropathy.** *Pain Manag.* 2015;5(4):285-296.
    <https://pubmed.ncbi.nlm.nih.gov/26087973/>

45. Robinson CR, Dougherty PM. **Spinal astrocyte gap junction and glutamate transporter expression contributes to a rat model of bortezomib-induced peripheral neuropathy.** *Neuroscience.* 2015;285:1-10.
    <https://pubmed.ncbi.nlm.nih.gov/25446351/>
    — reduced GLT-1 → an increase in `CENTS` gain.

46. Nashawi H, Masocha W, Edafiogho IO, Kombian SB. **Paclitaxel Causes Electrophysiological Changes in the Anterior Cingulate Cortex via Modulation of the γ-Aminobutyric Acid-ergic System.** *Med Princ Pract.* 2016;25(5):423-428.
    <https://pubmed.ncbi.nlm.nih.gov/27285838/>

47. Ossipov MH, Morimura K, Porreca F. **Descending pain modulation and chronification of pain.** *Curr Opin Support Palliat Care.* 2014;8(2):143-151.
    <https://pubmed.ncbi.nlm.nih.gov/24752199/>
    — descending noradrenergic inhibition — the site of action of duloxetine (`NATONE`).

---

## 10. Clinical pharmacokinetics — the model PK parameters

48. Graham MA, et al. **Clinical pharmacokinetics of oxaliplatin: a critical review.** *Clin Cancer Res.* 2000;6(4):1205-1218.
    <https://pubmed.ncbi.nlm.nih.gov/10778943/>
    — ultrafiltrable platinum Cmax ~0.81 µg/mL (85 mg/m², 2 h), CL 13.3 L/h, t½α 0.43 h, t½γ 391 h. The direct basis for the three-compartment PK of oxaliplatin.

49. Gamelin E, et al. **Cumulative pharmacokinetic study of oxaliplatin, administered every three weeks, combined with 5-fluorouracil in colorectal cancer patients.** *Clin Cancer Res.* 1997;3(6):891-899.
    <https://pubmed.ncbi.nlm.nih.gov/9815763/>
    — the long-term accumulation of tissue-bound platinum. The deep compartment (t½ ≈ 14 d).

50. Henningsson A, et al. **Mechanism-based pharmacokinetic model for paclitaxel.** *J Clin Oncol.* 2001;19(20):4065-4073.
    <https://pubmed.ncbi.nlm.nih.gov/11600609/>
    — paclitaxel three-compartment plus Michaelis-Menten elimination (Vmax ≈ 30 mg/h, Km ≈ 0.73 mg/L).

51. Joerger M, et al. **Population pharmacokinetics and pharmacodynamics of paclitaxel and carboplatin in ovarian cancer patients.** *Clin Cancer Res.* 2007;13(21):6410-6418.
    <https://pubmed.ncbi.nlm.nih.gov/17975154/>
    — Tc>0.05 µmol/L is the PD driving variable for toxicity. The `TCTHR` tracking variable.

52. Mielke S, et al. **Association of Paclitaxel Pharmacokinetics with the Development of Peripheral Neuropathy in Patients with Advanced Cancer.** *Clin Cancer Res.* 2005;11(13):4843-4850.
    <https://pubmed.ncbi.nlm.nih.gov/16000582/>
    — the association of time above threshold with neuropathy.

53. Moreau P, et al. **Subcutaneous versus intravenous administration of bortezomib in patients with relapsed multiple myeloma: a randomised, phase 3, non-inferiority study (MMY-3021).** *Lancet Oncol.* 2011;12(5):431-440.
    <https://pubmed.ncbi.nlm.nih.gov/21507715/>
    — **the key validation basis of the model.** SC against IV: Cmax 20.4 vs 223 ng/mL, AUC 155 vs 151 ng·h/mL (equivalent), grade ≥2 peripheral neuropathy 24% vs 41%, response rate identical.

54. Reece DE, et al. **Pharmacokinetic and pharmacodynamic study of two doses of bortezomib in patients with relapsed multiple myeloma.** *Cancer Chemother Pharmacol.* 2011;67(1):57-67.
    <https://pubmed.ncbi.nlm.nih.gov/20074322/>
    — the time course of whole-blood 20S proteasome inhibition (reversible, t½ ≈ 110 min).

55. Lantz RJ, et al. **Metabolism, excretion, and pharmacokinetics of duloxetine in healthy human subjects.** *Drug Metab Dispos.* 2003;31(9):1142-1150.
    <https://pubmed.ncbi.nlm.nih.gov/12920170/>
    — the PK of duloxetine (t½ ≈ 12 h, CYP1A2/2D6).

---

## 11. Pivotal trials — calibration and validation anchors

56. Grothey A, et al. (IDEA Collaboration) **Duration of Adjuvant Chemotherapy for Stage III Colon Cancer.** *N Engl J Med.* 2018;378(13):1177-1188.
    <https://pubmed.ncbi.nlm.nih.gov/29590544/>
    — **the primary calibration basis of the model.** FOLFOX for 3 months against 6 months, grade ≥2 neuropathy 16.6% vs 47.7%; 3-year DFS 83.1% vs 83.3% in the low-risk group and 62.7% vs 64.4% in the high-risk group. The neuropathy of the 3-month arm is an out-of-sample prediction target for the model.

57. André T, et al. **Oxaliplatin, fluorouracil, and leucovorin as adjuvant treatment for colon cancer (MOSAIC).** *N Engl J Med.* 2004;350(23):2343-2351.
    <https://pubmed.ncbi.nlm.nih.gov/15175436/>
    — grade 3 neuropathy 12.4% (during treatment).

58. André T, et al. **Improved overall survival with oxaliplatin, fluorouracil, and leucovorin as adjuvant treatment in stage II or III colon cancer in the MOSAIC trial.** *J Clin Oncol.* 2009;27(19):3109-3116.
    <https://pubmed.ncbi.nlm.nih.gov/19451431/>
    — **the basis for validating recovery.** Grade 3 neuropathy 12.4% → 1.1% at 12 months → 0.7% at 48 months.

59. Sparano JA, et al. (ECOG 1199) **Weekly paclitaxel in the adjuvant treatment of breast cancer.** *N Engl J Med.* 2008;358(16):1663-1671.
    <https://pubmed.ncbi.nlm.nih.gov/18420499/>
    — weekly paclitaxel 80 mg/m² against 175 mg/m² every three weeks: grade ≥2 sensory neuropathy 27% vs 20%. The q3w arm is an out-of-sample prediction target for the model.

60. de Gramont A, et al. **OPTIMOX1: a randomized study of FOLFOX4 or FOLFOX7 with oxaliplatin in a stop-and-Go fashion in advanced colorectal cancer.** *J Clin Oncol.* 2006;24(3):394-400.
    <https://pubmed.ncbi.nlm.nih.gov/16421419/>
    — the reduction in neurotoxicity with a stop-and-go strategy. Scenario L of the model.

61. Smith EML, et al. **Effect of duloxetine on pain, function, and quality of life among patients with chemotherapy-induced painful peripheral neuropathy: a randomized clinical trial.** *JAMA.* 2013;309(13):1359-1367.
    <https://pubmed.ncbi.nlm.nih.gov/23549581/>
    — **the calibration basis for duloxetine.** A fall in average pain over 5 weeks of 1.06 against 0.34 on placebo (BPI 0-10).

62. Rao RD, et al. **Efficacy of gabapentin in the management of chemotherapy-induced peripheral neuropathy: a phase 3 randomized, double-blind, placebo-controlled, crossover trial (N00C3).** *Cancer.* 2007;110(9):2110-2118.
    <https://pubmed.ncbi.nlm.nih.gov/17853395/>
    — the negative gabapentin result.

63. Hershman DL, et al. **Randomized double-blind placebo-controlled trial of acetyl-L-carnitine for the prevention of taxane-induced neuropathy in women undergoing adjuvant breast cancer therapy (SWOG S0715).** *J Clin Oncol.* 2013;31(20):2627-2633.
    <https://pubmed.ncbi.nlm.nih.gov/23733756/>
    — acetyl-L-carnitine made things **worse**. The harmful node of cluster 13.

64. Loprinzi CL, et al. **Phase III randomized, placebo-controlled, double-blind study of intravenous calcium and magnesium to prevent oxaliplatin-induced sensory neurotoxicity (N08CB/Alliance).** *J Clin Oncol.* 2014;32(10):997-1005.
    <https://pubmed.ncbi.nlm.nih.gov/24297951/>
    — the negative Ca/Mg result — the limits of the oxalate chelation hypothesis.

65. Hanai A, et al. **Effects of Cryotherapy on Objective and Subjective Symptoms of Paclitaxel-Induced Neuropathy: Prospective Self-Controlled Trial.** *J Natl Cancer Inst.* 2018;110(2):141-148.
    <https://pubmed.ncbi.nlm.nih.gov/29924336/>
    — the preventive effect of cryotherapy. The model's `CRYO` (the fraction of delivery blocked).

66. Michel P, et al. / Sundar R, et al. **Limb Hypothermia for Preventing Paclitaxel-Induced Peripheral Neuropathy in Breast Cancer Patients: A Pilot Study.** *Front Oncol.* 2017;6:274.
    <https://pubmed.ncbi.nlm.nih.gov/28123993/>

67. Kanbayashi Y, et al. **Effects of compression therapy for prevention of taxane-induced peripheral neuropathy.** — Tsuyuki S, et al. **Evaluation of the effect of compression therapy using surgical gloves on nanoparticle albumin-bound paclitaxel-induced peripheral neuropathy: a phase II multicenter study.** *Breast Cancer Res Treat.* 2016;160(1):61-67.
    <https://pubmed.ncbi.nlm.nih.gov/27620884/>

68. Kleckner IR, et al. **Effects of exercise during chemotherapy on chemotherapy-induced peripheral neuropathy: a multicenter, randomized controlled trial.** *Support Care Cancer.* 2018;26(4):1019-1028.
    <https://pubmed.ncbi.nlm.nih.gov/29028095/>

---

## 12. Risk factors and pharmacogenomics

69. Argyriou AA, et al. **Diabetes mellitus as a risk factor for chemotherapy-induced peripheral neuropathy.** — Uwah AN, et al. **The effect of diabetes on oxaliplatin-induced peripheral neuropathy.** *Clin Colorectal Cancer.* 2012;11(4):275-279.
    <https://pubmed.ncbi.nlm.nih.gov/22682697/>
    — the reduced baseline axonal reserve of the diabetic patient. The model's `AXON0` and `RESERVE` covariates.

70. Baldwin RM, et al. **A genome-wide association study identifies novel loci for paclitaxel-induced sensory peripheral neuropathy in CALGB 40101.** *Clin Cancer Res.* 2012;18(18):5099-5109.
    <https://pubmed.ncbi.nlm.nih.gov/22843789/>
    — EPHA5 and others. Cluster 12.

71. Hertz DL, et al. **CYP2C8*3 increases risk of neuropathy in breast cancer patients treated with paclitaxel.** *Ann Oncol.* 2013;24(6):1472-1478.
    <https://pubmed.ncbi.nlm.nih.gov/23413280/>

72. Chua KC, et al. **Genetic variants of CYP3A4 and CIPN risk.** — Marcath LA, et al. **Comprehensive assessment of cytochromes P450 and transporter genetics with endocrine treatment and taxane neuropathy.** *Pharmacogenet Genomics.* 2019;29(4):79-86.
    <https://pubmed.ncbi.nlm.nih.gov/30920425/>

73. Chaudhry V, Chaudhry M, Crawford TO, Simmons-O'Brien E, Griffin JW. **Toxic neuropathy in patients with pre-existing neuropathy.** *Neurology.* 2003;60(2):337-340.
    <https://pubmed.ncbi.nlm.nih.gov/12552056/>
    — the vulnerability of those who already carry a neuropathy such as CMT. The `CMT1A carrier` scenario.

74. Lee JJ, Swain SM. **Peripheral neuropathy induced by microtubule-stabilizing agents.** *J Clin Oncol.* 2006;24(10):1633-1642.
    <https://pubmed.ncbi.nlm.nih.gov/16575015/>

---

## 13. Biomarkers and objective measures

75. Karteri S, et al. **Prospectively Assessing Serum Neurofilament Light Chain Levels As A Biomarker Of Paclitaxel-Induced Peripheral Neurotoxicity In Breast Cancer Patients.** *J Peripher Nerv Syst.* 2022;27(1):166-174.
    <https://pubmed.ncbi.nlm.nih.gov/35112433/>
    — **the direct basis for the `NFL` compartment.** The rise in NfL precedes the change in the clinical scales.

76. Huehnchen P, et al. **Neurofilament proteins as a potential biomarker in chemotherapy-induced polyneuropathy.** *JCI Insight.* 2022;7(6):e154395.
    <https://pubmed.ncbi.nlm.nih.gov/35133982/>
    — NfL predicts CIPN early. The basis for the lead-time analysis of the model.

77. Meregalli C, et al. **Neurofilament light chain as disease biomarker in a rodent model of chemotherapy induced peripheral neuropathy.** *Exp Neurol.* 2018;307:129-132.
    <https://pubmed.ncbi.nlm.nih.gov/29908886/>

78. Lauria G, et al. **European Federation of Neurological Societies/Peripheral Nerve Society Guideline on the use of skin biopsy in the diagnosis of small fiber neuropathy.** *Eur J Neurol.* 2010;17(7):903-912.
    <https://pubmed.ncbi.nlm.nih.gov/20642627/>
    — the normal range of IENFD (the model's reference value of 7.0 fibres/mm).

79. Postma TJ, et al. **The development of an EORTC quality of life questionnaire to assess chemotherapy-induced peripheral neuropathy: the QLQ-CIPN20.** *Eur J Cancer.* 2005;41(8):1135-1139.
    <https://pubmed.ncbi.nlm.nih.gov/15911236/>
    — the definition of the CIPN20 scale.

80. Cavaletti G, et al. **The Total Neuropathy Score as an assessment tool for grading the course of chemotherapy-induced peripheral neurotoxicity: comparison with the National Cancer Institute-Common Toxicity Scale.** *J Peripher Nerv Syst.* 2007;12(3):210-215.
    <https://pubmed.ncbi.nlm.nih.gov/17868248/>
    — the relationship of TNSc to the CTCAE grade — the mapping of the model's grade thresholds.

81. Alberti P, et al. **Physician-assessed and patient-reported outcome measures in chemotherapy-induced sensory peripheral neurotoxicity: two sides of the same coin.** *Ann Oncol.* 2014;25(1):257-264.
    <https://pubmed.ncbi.nlm.nih.gov/24356636/>
    — the discordance between the clinician's grade and the patient's report — why the model separates CTCAE (structure-weighted) from CIPN20 (which includes symptoms).

---

## 14. QSP and modelling methodology

82. Elmeliegy M, et al. **Towards better understanding of CIPN: mechanistic modelling.** — Ballesta A, et al. **A mixed-effects model of the dynamics of chemotherapy-induced peripheral neuropathy.** *CPT Pharmacometrics Syst Pharmacol.* 2019;. (methodology comparator)
    <https://pubmed.ncbi.nlm.nih.gov/31215775/>

83. Baaz M, Cardilin T, Jirstrand M. **Model-Based Prediction of Progression-Free Survival for Combination Therapies in Oncology.** *CPT Pharmacometrics Syst Pharmacol.* 2023;12(9):1227-1237.
    <https://pubmed.ncbi.nlm.nih.gov/37448175/>
    — the conversion of tumour burden into a survival endpoint (the model's Poisson zero-clonogen DFS mapping).

84. Norton L, Simon R. **Tumor size, sensitivity to therapy, and design of treatment schedules.** *Cancer Treat Rep.* 1977;61(7):1307-1317.
    <https://pubmed.ncbi.nlm.nih.gov/589597/>
    — log-kill theory — the basis for the tumour axis.

85. Baker SD, et al. **Relationship of systemic exposure to unbound paclitaxel and neutropenia.** *Clin Pharmacol Ther.* 2006;80(1):72-83.
    <https://pubmed.ncbi.nlm.nih.gov/16815320/>

86. Elharrar X, et al. / Cavaletti G, Marmiroli P. **Management of oxaliplatin-induced peripheral sensory neuropathy.** *Cancers (Basel).* 2020;12(6):1370.
    <https://pubmed.ncbi.nlm.nih.gov/32471249/>

87. Gewandter JS, et al. **Trial designs for chemotherapy-induced peripheral neuropathy prevention: ACTTION recommendations.** *Neurology.* 2018;91(9):403-413.
    <https://pubmed.ncbi.nlm.nih.gov/30054438/>
    — the design of prevention trials — the basis for the comparative scenario design of the model.

88. Kolb NA, et al. **The Association of Chemotherapy-Induced Peripheral Neuropathy Symptoms and the Risk of Falling.** *JAMA Neurol.* 2016;73(7):860-866.
    <https://pubmed.ncbi.nlm.nih.gov/27183099/>
    — the fall risk endpoint (cluster 15).

89. Winters-Stone KM, et al. **Falls, Functioning, and Disability Among Women With Persistent Symptoms of Chemotherapy-Induced Peripheral Neuropathy.** *J Clin Oncol.* 2017;35(23):2604-2612.
    <https://pubmed.ncbi.nlm.nih.gov/28586243/>

90. Timmins HC, et al. **Taxane-induced peripheral neuropathy: differences in patient report and objective assessment.** *Support Care Cancer.* 2020;28(9):4459-4466.
    <https://pubmed.ncbi.nlm.nih.gov/31925533/>

---

## 15. Therapeutic strategies and future directions

91. Kim JH, et al. **HDAC6 inhibitors rescued the defective axonal mitochondrial movement in motor neurons.** — Krukowski K, et al. **HDAC6 inhibition effectively reverses chemotherapy-induced peripheral neuropathy.** *Pain.* 2017;158(6):1126-1137.
    <https://pubmed.ncbi.nlm.nih.gov/28267067/>
    — restoration of axonal transport by HDAC6 inhibition — an investigational node of cluster 13.

92. Chen Y, et al. **S1PR1 antagonist as a novel target for CIPN.** — Stockstill K, et al. **Dysregulation of sphingolipid metabolism contributes to bortezomib-induced neuropathic pain.** *J Exp Med.* 2018;215(5):1301-1313.
    <https://pubmed.ncbi.nlm.nih.gov/29703731/>

93. Smith EML, et al. **Patient-reported outcome measures in CIPN trials.** — Hershman DL, et al. **Comparison of physician- and patient-reported symptoms of neuropathy.** *Breast Cancer Res Treat.* 2011;125(3):767-774.
    <https://pubmed.ncbi.nlm.nih.gov/21128110/>

94. Pachman DR, et al. **Scrambler therapy for chemotherapy neuropathy: a randomized phase II pilot trial.** *Support Care Cancer.* 2019;27(11):4123-4131.
    <https://pubmed.ncbi.nlm.nih.gov/30822010/>

95. Molassiotis A, et al. **Acupuncture for chemotherapy-induced peripheral neuropathy: a pilot randomised controlled trial.** *Clin Trials Regul Sci Cancer.* 2019. / Lu W, et al. **Acupuncture for Chemotherapy-Induced Peripheral Neuropathy in Breast Cancer Survivors: A Randomized Controlled Pilot Trial.** *Oncologist.* 2020;25(4):310-318.
    <https://pubmed.ncbi.nlm.nih.gov/32297442/>

96. Dorsey SG, et al. **NCI Clinical Trials Planning Meeting for prevention and treatment of chemotherapy-induced peripheral neuropathy.** *J Natl Cancer Inst.* 2019;111(6):531-537.
    <https://pubmed.ncbi.nlm.nih.gov/30715378/>

97. Cavaletti G, Marmiroli P. **Chemotherapy-induced peripheral neurotoxicity: A multifaceted, still unsolved issue.** *J Peripher Nerv Syst.* 2019;24 Suppl 2:S6-S12.
    <https://pubmed.ncbi.nlm.nih.gov/31647153/>

98. Tofthagen C, et al. **Falls in persons with chemotherapy-induced peripheral neuropathy.** *Support Care Cancer.* 2012;20(3):583-589.
    <https://pubmed.ncbi.nlm.nih.gov/21380613/>

99. Beijers AJM, et al. **Chemotherapy-induced peripheral neuropathy and impact on quality of life 6 months after treatment with chemotherapy.** *J Community Support Oncol.* 2014;12(11):401-406.
    <https://pubmed.ncbi.nlm.nih.gov/25844425/>

100. Grothey A. **Clinical management of oxaliplatin-associated neurotoxicity.** *Clin Colorectal Cancer.* 2005;5 Suppl 1:S38-46.
     <https://pubmed.ncbi.nlm.nih.gov/15871765/>
     — the clinical logic of dose modification strategies — the background to the therapeutic index analysis of §11 of the model.

---

## 16. Tools

101. Elmokadem A, Riggs MM, Baron KT. **Quantitative Systems Pharmacology and Physiologically-Based Pharmacokinetic Modeling With mrgsolve: A Hands-On Tutorial.** *CPT Pharmacometrics Syst Pharmacol.* 2019;8(12):883-893.
     <https://pubmed.ncbi.nlm.nih.gov/31654500/>

102. Baron KT. **mrgsolve: Simulate from ODE-Based Models.** <https://mrgsolve.org/>

---

**102 references in total** (with PubMed links). The provenance of the model parameters is also given alongside in the
`$PARAM` block comments of `cipn_mrgsolve_model.R` and in the calibration section of `cipn_reference_model.py`.
