# Periodontitis QSP Model — References

Each entry notes which part of the model it supports.
Each entry notes which part of the model it supports.

---

## 1. The framework of the aetiology — the polymicrobial synergy and dysbiosis (PSD) model

The basis for the microbial compartment structure of this model (a three-layer structure separating the commensals, the dysbiotic organisms, and the keystone pathogen) and
for the low-concentration, high-leverage term of the form `KEY = 1 + EKEY·PG/(KMKEY+PG)`.

1. Hajishengallis G, Darveau RP, Curtis MA. **The keystone-pathogen hypothesis.** Nat Rev Microbiol. 2012;10(10):717-25. <https://pubmed.ncbi.nlm.nih.gov/22941505/>
2. Hajishengallis G, Lamont RJ. **Beyond the red complex and into more complexity: the polymicrobial synergy and dysbiosis (PSD) model of periodontal disease etiology.** Mol Oral Microbiol. 2012;27(6):409-19. <https://pubmed.ncbi.nlm.nih.gov/23134607/>
3. Hajishengallis G, Liang S, Payne MA, et al. **Low-abundance biofilm species orchestrates inflammatory periodontal disease through the commensal microbiota and complement.** Cell Host Microbe. 2011;10(5):497-506. <https://pubmed.ncbi.nlm.nih.gov/22036469/>
4. Socransky SS, Haffajee AD, Cugini MA, Smith C, Kent RL. **Microbial complexes in subgingival plaque.** J Clin Periodontol. 1998;25(2):134-44. <https://pubmed.ncbi.nlm.nih.gov/9495612/>
5. Darveau RP. **Periodontitis: a polymicrobial disruption of host homeostasis.** Nat Rev Microbiol. 2010;8(7):481-90. <https://pubmed.ncbi.nlm.nih.gov/20514045/>
6. Lamont RJ, Koo H, Hajishengallis G. **The oral microbiota: dynamic communities and host interactions.** Nat Rev Microbiol. 2018;16(12):745-59. <https://pubmed.ncbi.nlm.nih.gov/30301974/>
7. Abusleme L, Dupuy AK, Dutzan N, et al. **The subgingival microbiome in health and periodontitis and its relationship with community biomass and inflammation.** ISME J. 2013;7(5):1016-25. <https://pubmed.ncbi.nlm.nih.gov/23303375/>

## 2. The inflammophilic nutritional feedback — the engine of the model

The basis for the `NUTR` (GCF-derived nutrient) → dysbiotic organism growth term. That this feedback
makes the loop gain > 1 and so gives rise to bistability is the first organising principle of this model.

8. Hajishengallis G. **The inflammophilic character of the periodontitis-associated microbiota.** Mol Oral Microbiol. 2014;29(6):248-57. <https://pubmed.ncbi.nlm.nih.gov/24976068/>
9. Herrero ER, Fernandes S, Verspecht T, et al. **Dysbiotic biofilms deregulate the periodontal inflammatory response.** J Dent Res. 2018;97(5):547-55. <https://pubmed.ncbi.nlm.nih.gov/29303019/>
10. Smalley JW, Olczak T. **Haem acquisition mechanisms of Porphyromonas gingivalis — strategies used in a polymicrobial community in a haem-limited host environment.** Mol Oral Microbiol. 2017;32(1):1-23. <https://pubmed.ncbi.nlm.nih.gov/26662717/>
11. Griffiths GS. **Formation, collection and significance of gingival crevice fluid.** Periodontol 2000. 2003;31:32-42. <https://pubmed.ncbi.nlm.nih.gov/12656994/>

## 3. Keystone pathogen virulence — gingipains, LPS heterogeneity, tissue invasion

The basis for `GING` (gingipain activity), the C5-convertase-like activity of the gingipains (C3-independent C5a generation), and the
tissue invasion fraction (`TISS_FRAC`).

12. Potempa J, Sroka A, Imamura T, Travis J. **Gingipains, the major cysteine proteinases and virulence factors of Porphyromonas gingivalis.** Curr Protein Pept Sci. 2003;4(6):397-407. <https://pubmed.ncbi.nlm.nih.gov/14683426/>
13. Wingrove JA, DiScipio RG, Chen Z, Potempa J, Travis J, Hugli TE. **Activation of complement components C3 and C5 by a cysteine proteinase (gingipain-1) from Porphyromonas gingivalis.** J Biol Chem. 1992;267(26):18902-7. <https://pubmed.ncbi.nlm.nih.gov/1527018/>
14. Coats SR, Jones JW, Do CT, et al. **Human Toll-like receptor 4 responses to P. gingivalis are regulated by lipid A 1- and 4'-phosphatase activities.** Cell Microbiol. 2009;11(11):1587-99. <https://pubmed.ncbi.nlm.nih.gov/19552698/>
15. Lamont RJ, Chan A, Belton CM, Izutsu KT, Vasel D, Weinberg A. **Porphyromonas gingivalis invasion of gingival epithelial cells.** Infect Immun. 1995;63(10):3878-85. <https://pubmed.ncbi.nlm.nih.gov/7558295/>
16. Wegner N, Wait R, Sroka A, et al. **Peptidylarginine deiminase from Porphyromonas gingivalis citrullinates human fibrinogen and alpha-enolase.** Arthritis Rheum. 2010;62(9):2662-72. <https://pubmed.ncbi.nlm.nih.gov/20506214/>
17. Kachlany SC. **Aggregatibacter actinomycetemcomitans leukotoxin: from threat to therapy.** J Dent Res. 2010;89(6):561-70. <https://pubmed.ncbi.nlm.nih.gov/20200418/>
18. Haubek D, Ennibi OK, Poulsen K, Væth M, Poulsen S, Kilian M. **Risk of aggressive periodontitis in adolescent carriers of the JP2 clone of A. actinomycetemcomitans in Morocco: a prospective longitudinal cohort study.** Lancet. 2008;371(9608):237-42. <https://pubmed.ncbi.nlm.nih.gov/18207019/>

## 4. Complement–TLR2 crosstalk and the subversion of neutrophil killing — the second organising principle of the model

The basis for holding killing capacity (`PKILL`) as a **state variable** rather than a constant, for the Hill form of the `SUBV` term,
and for the structure making bacterial clearance proportional to `PMN × PKILL`.

19. Maekawa T, Krauss JL, Abe T, et al. **Porphyromonas gingivalis manipulates complement and TLR signaling to uncouple bacterial clearance from inflammation and promote dysbiosis.** Cell Host Microbe. 2014;15(6):768-78. <https://pubmed.ncbi.nlm.nih.gov/24922578/>
20. Wang M, Krauss JL, Domon H, et al. **Microbial hijacking of complement-Toll-like receptor crosstalk.** Sci Signal. 2010;3(109):ra11. <https://pubmed.ncbi.nlm.nih.gov/20159852/>
21. Hajishengallis G, Lambris JD. **Microbial manipulation of receptor crosstalk in innate immunity.** Nat Rev Immunol. 2011;11(3):187-200. <https://pubmed.ncbi.nlm.nih.gov/21331082/>
22. Hajishengallis G, Lambris JD. **Complement and dysbiosis in periodontal disease.** Immunobiology. 2012;217(11):1111-6. <https://pubmed.ncbi.nlm.nih.gov/22964237/>
23. Hajishengallis G, Chavakis T, Hajishengallis E, Lambris JD. **Neutrophil homeostasis and inflammation: novel paradigms from studying periodontitis.** J Leukoc Biol. 2015;98(4):539-48. <https://pubmed.ncbi.nlm.nih.gov/25548253/>
24. Hajishengallis G, Moutsopoulos NM, Hajishengallis E, Chavakis T. **Immune and regulatory functions of neutrophils in inflammatory bone loss.** Semin Immunol. 2016;28(2):146-58. <https://pubmed.ncbi.nlm.nih.gov/26936034/>

## 5. The paradox of neutrophil deficiency — LAD-I and IL-17

The basis for the structure in which the `LADI` switch abolishes neutrophil recruitment and yet worsens bone loss through IL-17.

25. Moutsopoulos NM, Konkel J, Sarmadi M, et al. **Defective neutrophil recruitment in leukocyte adhesion deficiency type I disease causes local IL-17-driven inflammatory bone loss.** Sci Transl Med. 2014;6(229):229ra40. <https://pubmed.ncbi.nlm.nih.gov/24670684/>
26. Moutsopoulos NM, Zerbe CS, Wild T, et al. **Interleukin-12 and interleukin-23 blockade in leukocyte adhesion deficiency type 1.** N Engl J Med. 2017;376(12):1141-6. <https://pubmed.ncbi.nlm.nih.gov/28328334/>
27. Dutzan N, Kajikawa T, Abusleme L, et al. **A dysbiotic microbiome triggers TH17 cells to mediate oral mucosal immunopathology in mice and humans.** Sci Transl Med. 2018;10(463):eaat0797. <https://pubmed.ncbi.nlm.nih.gov/30333238/>
28. Dutzan N, Abusleme L, Bridgeman H, et al. **On-going mechanical damage from mastication drives homeostatic Th17 cell responses at the oral barrier.** Immunity. 2017;46(1):133-47. <https://pubmed.ncbi.nlm.nih.gov/28087239/>

## 6. RANKL/OPG and bone loss — the ratchet of the model

The basis for the cell-source weights of `RANKL_T` (Th17, B and plasma cells, stromal cells), the `RATIO^HOC` osteoclast drive,
and the structure that makes `BLOSS` a monotonically increasing integral.

29. Boyle WJ, Simonet WS, Lacey DL. **Osteoclast differentiation and activation.** Nature. 2003;423(6937):337-42. <https://pubmed.ncbi.nlm.nih.gov/12748652/>
30. Kawai T, Matsuyama T, Hosokawa Y, et al. **B and T lymphocytes are the primary sources of RANKL in the bone resorptive lesion of periodontal disease.** Am J Pathol. 2006;169(3):987-98. <https://pubmed.ncbi.nlm.nih.gov/16936272/>
31. Belibasakis GN, Bostanci N. **The RANKL-OPG system in clinical periodontology.** J Clin Periodontol. 2012;39(3):239-48. <https://pubmed.ncbi.nlm.nih.gov/22092994/>
32. Graves DT, Cochran D. **The contribution of interleukin-1 and tumor necrosis factor to periodontal tissue destruction.** J Periodontol. 2003;74(3):391-401. <https://pubmed.ncbi.nlm.nih.gov/12710761/>
33. Waerhaug J. **The angular bone defect and its relationship to trauma from occlusion and downgrowth of subgingival plaque.** J Clin Periodontol. 1979;6(2):61-82. <https://pubmed.ncbi.nlm.nih.gov/286657/>
34. Page RC, Schroeder HE. **Pathogenesis of inflammatory periodontal disease. A summary of current work.** Lab Invest. 1976;34(3):235-49. <https://pubmed.ncbi.nlm.nih.gov/765622/>

## 7. MMP–TIMP imbalance and connective tissue destruction

The basis for computing `MMP8A` (active MMP-8) as **activity** rather than as a total amount, and for having doxycycline inhibit
the activity rather than the synthesis.

35. Sorsa T, Tjäderhane L, Konttinen YT, et al. **Matrix metalloproteinases: contribution to pathogenesis, diagnosis and treatment of periodontal inflammation.** Ann Med. 2006;38(5):306-21. <https://pubmed.ncbi.nlm.nih.gov/16938801/>
36. Golub LM, Lee HM, Ryan ME, Giannobile WV, Payne J, Sorsa T. **Tetracyclines inhibit connective tissue breakdown by multiple non-antimicrobial mechanisms.** Adv Dent Res. 1998;12(2):12-26. <https://pubmed.ncbi.nlm.nih.gov/9972117/>
37. Sorsa T, Alassiri S, Grigoriadis A, et al. **Active MMP-8 (aMMP-8) as a grading and staging biomarker in the periodontitis classification.** Diagnostics (Basel). 2020;10(2):61. <https://pubmed.ncbi.nlm.nih.gov/31979091/>
38. Uitto VJ, Overall CM, McCulloch C. **Proteolytic host cell enzymes in gingival crevice fluid.** Periodontol 2000. 2003;31:77-104. <https://pubmed.ncbi.nlm.nih.gov/12656997/>

## 8. Failure of resolution and the SPMs

The basis for the `RVELOC` (resolvin E1) compartment, and for the **non-immunosuppressive** mode of action that combines
suppression of neutrophil recruitment with enhanced efficient clearance of apoptotic cells (efferocytosis).

39. Serhan CN. **Pro-resolving lipid mediators are leads for resolution physiology.** Nature. 2014;510(7503):92-101. <https://pubmed.ncbi.nlm.nih.gov/24899309/>
40. Hasturk H, Kantarci A, Goguet-Surmenian E, et al. **Resolvin E1 regulates inflammation at the cellular and tissue level and restores tissue homeostasis in vivo.** J Immunol. 2007;179(10):7021-9. <https://pubmed.ncbi.nlm.nih.gov/17982093/>
41. Van Dyke TE. **Pro-resolving mediators in the regulation of periodontal disease.** Mol Aspects Med. 2017;58:21-36. <https://pubmed.ncbi.nlm.nih.gov/28483532/>
42. Van Dyke TE, Hasturk H, Kantarci A, et al. **Proresolving nanomedicines activate bone regeneration in periodontitis.** J Dent Res. 2015;94(1):148-56. <https://pubmed.ncbi.nlm.nih.gov/25389003/>

## 9. Mechanical therapy — the evidence that the depth dependence is anatomical

The direct data source for the `SRP_EFF` logistic (`PPD50_SRP = 5.8 mm`, `PPD_W = 1.5 mm`).
The model predicts 82% / 63% / 34% calculus removal at 3.0 / 4.5 / 6.5 mm respectively, and
Waerhaug's measurements are 83% / 61% / 32%.

43. Waerhaug J. **Healing of the dento-epithelial junction following subgingival plaque control. II. As observed on extracted teeth.** J Periodontol. 1978;49(3):119-34. <https://pubmed.ncbi.nlm.nih.gov/347090/>
44. Rabbani GM, Ash MM, Caffesse RG. **The effectiveness of subgingival scaling and root planing in calculus removal.** J Periodontol. 1981;52(3):119-23. <https://pubmed.ncbi.nlm.nih.gov/6939862/>
45. Stambaugh RV, Dragoo M, Smith DM, Carasali L. **The limits of subgingival scaling.** Int J Periodontics Restorative Dent. 1981;1(5):30-41. <https://pubmed.ncbi.nlm.nih.gov/7047430/>
46. Cobb CM. **Non-surgical pocket therapy: mechanical.** Ann Periodontol. 1996;1(1):443-90. <https://pubmed.ncbi.nlm.nih.gov/9118261/>
47. Cobb CM. **Clinical significance of non-surgical periodontal therapy: an evidence-based perspective of scaling and root planing.** J Clin Periodontol. 2002;29 Suppl 2:6-16. <https://pubmed.ncbi.nlm.nih.gov/12010523/>

## 10. The critical probing depth — an emergent phenomenon the model has to reproduce

The model did not use these values as an objective function; from the single trade-off "surgery does not change the host parameters, only the accessibility
(`PPD50`), and pays a fixed cost in recession and attachment loss" it produces
2.13 mm for SRP and 3.43 mm for surgery (the literature values being 2.9 / 4.2 mm).

48. Lindhe J, Socransky SS, Nyman S, Haffajee A, Westfelt E. **"Critical probing depths" in periodontal therapy.** J Clin Periodontol. 1982;9(4):323-36. <https://pubmed.ncbi.nlm.nih.gov/6957467/>
49. Heitz-Mayfield LJ, Trombelli L, Heitz F, Needleman I, Moles D. **A systematic review of the effect of surgical debridement vs. non-surgical debridement for the treatment of chronic periodontitis.** J Clin Periodontol. 2002;29 Suppl 3:92-102. <https://pubmed.ncbi.nlm.nih.gov/12787211/>
50. Badersten A, Nilvéus R, Egelberg J. **Effect of nonsurgical periodontal therapy. II. Severely advanced periodontitis.** J Clin Periodontol. 1984;11(1):63-76. <https://pubmed.ncbi.nlm.nih.gov/6363463/>
51. Sanz M, Herrera D, Kebschull M, et al. **Treatment of stage I-III periodontitis — the EFP S3 level clinical practice guideline.** J Clin Periodontol. 2020;47 Suppl 22:4-60. <https://pubmed.ncbi.nlm.nih.gov/32383274/>

## 11. Adjunctive antibiotics

The basis for `EMAX_AMX`/`EMAX_MTZ`/`SYN_AM` and for the model's prediction that "the benefit concentrates in the deep pockets".

52. van Winkelhoff AJ, Rodenburg JP, Goené RJ, Abbas F, Winkel EG, de Graaff J. **Metronidazole plus amoxycillin in the treatment of Actinobacillus actinomycetemcomitans associated periodontitis.** J Clin Periodontol. 1989;16(2):128-31. <https://pubmed.ncbi.nlm.nih.gov/2921374/>
53. Feres M, Soares GM, Mendes JA, et al. **Metronidazole alone or with amoxicillin as adjuncts to non-surgical treatment of chronic periodontitis: a 1-year double-blinded, placebo-controlled, randomized clinical trial.** J Clin Periodontol. 2012;39(12):1149-58. <https://pubmed.ncbi.nlm.nih.gov/23016867/>
54. Harks I, Koch R, Eickholz P, et al. **Is progression of periodontitis relevantly influenced by systemic antibiotics? A clinical randomized trial.** J Clin Periodontol. 2015;42(9):832-42. <https://pubmed.ncbi.nlm.nih.gov/26250060/>
55. Teughels W, Feres M, Oud V, Martín C, Matesanz P, Herrera D. **Adjunctive effect of systemic antimicrobials in periodontitis therapy: a systematic review and meta-analysis.** J Clin Periodontol. 2020;47 Suppl 22:257-81. <https://pubmed.ncbi.nlm.nih.gov/31994207/>
56. Herrera D, Matesanz P, Martín C, Oud V, Feres M, Teughels W. **Adjunctive effect of locally delivered antimicrobials in periodontitis therapy: a systematic review and meta-analysis.** J Clin Periodontol. 2020;47 Suppl 22:239-56. <https://pubmed.ncbi.nlm.nih.gov/31912531/>
57. Williams RC, Paquette DW, Offenbacher S, et al. **Treatment of periodontitis by local administration of minocycline microspheres: a controlled trial.** J Periodontol. 2001;72(11):1535-44. <https://pubmed.ncbi.nlm.nih.gov/11759865/>
58. Jepsen K, Jepsen S. **Antibiotics/antimicrobials: systemic and local administration in the therapy of mild to moderately advanced periodontitis.** Periodontol 2000. 2016;71(1):82-112. <https://pubmed.ncbi.nlm.nih.gov/27045432/>

## 12. Host modulation therapy

The basis for `SDD` (sub-antimicrobial dose doxycycline), the C3 inhibitor AMY-101, and the anti-RANKL agents.
The model predicts that the maximum effect of AMY-101 is structurally limited by the gingipain-derived C5a.

59. Caton JG, Ciancio SG, Blieden TM, et al. **Treatment with subantimicrobial dose doxycycline improves the efficacy of scaling and root planing in patients with adult periodontitis.** J Periodontol. 2000;71(4):521-32. <https://pubmed.ncbi.nlm.nih.gov/10807113/>
60. Preshaw PM, Hefti AF, Novak MJ, et al. **Subantimicrobial dose doxycycline enhances the efficacy of scaling and root planing in chronic periodontitis: a multicenter trial.** J Periodontol. 2004;75(8):1068-76. <https://pubmed.ncbi.nlm.nih.gov/15455734/>
61. Golub LM, Lee HM. **Periodontal therapeutics: current host-modulation agents and future directions.** Periodontol 2000. 2020;82(1):186-204. <https://pubmed.ncbi.nlm.nih.gov/31850625/>
62. Hasturk H, Hajishengallis G, Lambris JD, Mastellos DC, Yancopoulou D. **Phase IIa clinical trial of complement C3 inhibitor AMY-101 in adults with periodontal inflammation.** J Clin Invest. 2021;131(23):e152973. <https://pubmed.ncbi.nlm.nih.gov/34623331/>
63. Maekawa T, Abe T, Hajishengallis E, et al. **Genetic and intervention studies implicating complement C3 as a major target for the treatment of periodontitis.** J Immunol. 2014;192(12):6020-7. <https://pubmed.ncbi.nlm.nih.gov/24808362/>
64. Hajishengallis G, Kajikawa T, Hajishengallis E, et al. **Complement-dependent mechanisms and interventions in periodontal disease.** Front Immunol. 2019;10:406. <https://pubmed.ncbi.nlm.nih.gov/30915073/>
65. Reid IR, Billington EO. **Drug therapy for osteoporosis in older adults.** Lancet. 2022;399(10329):1080-92. <https://pubmed.ncbi.nlm.nih.gov/35279261/>
66. Ruggiero SL, Dodson TB, Aghaloo T, Carlson ER, Ward BB, Kademani D. **AAOMS position paper on medication-related osteonecrosis of the jaws — 2022 update.** J Oral Maxillofac Surg. 2022;80(5):920-43. <https://pubmed.ncbi.nlm.nih.gov/35300956/>

## 13. Risk-modifying factors such as smoking and diabetes

The basis for the `SMOKE` (impaired neutrophil function plus masking of BOP) and `HG_GAIN` (impairment of killing by hyperglycaemia) terms.

67. Tomar SL, Asma S. **Smoking-attributable periodontitis in the United States: findings from NHANES III.** J Periodontol. 2000;71(5):743-51. <https://pubmed.ncbi.nlm.nih.gov/10872955/>
68. Dietrich T, Bernimoulin JP, Glynn RJ. **The effect of cigarette smoking on gingival bleeding.** J Periodontol. 2004;75(1):16-22. <https://pubmed.ncbi.nlm.nih.gov/15025213/>
69. Labriola A, Needleman I, Moles DR. **Systematic review of the effect of smoking on nonsurgical periodontal therapy.** Periodontol 2000. 2005;37:124-37. <https://pubmed.ncbi.nlm.nih.gov/15655029/>
70. Lalla E, Papapanou PN. **Diabetes mellitus and periodontitis: a tale of two common interrelated diseases.** Nat Rev Endocrinol. 2011;7(12):738-48. <https://pubmed.ncbi.nlm.nih.gov/21709707/>
71. Graves DT, Ding Z, Yang Y. **The impact of diabetes on periodontal diseases.** Periodontol 2000. 2020;82(1):214-24. <https://pubmed.ncbi.nlm.nih.gov/31850631/>
72. Löe H, Ånerud Å, Boysen H, Morrison E. **Natural history of periodontal disease in man. Rapid, moderate and no loss of attachment in Sri Lankan labourers 14 to 46 years of age.** J Clin Periodontol. 1986;13(5):431-45. <https://pubmed.ncbi.nlm.nih.gov/3487557/>

## 14. Systemic spread — PISA, HbA1c, endothelial function

The basis for the `PISA` formula, the `IL6S → CRP → INSR → HBA1C` chain, and the non-monotonic behaviour in which
FMD **transiently worsens** 24 hours after intensive treatment and then improves at 6 months.

73. Nesse W, Abbas F, van der Ploeg I, Spijkervet FK, Dijkstra PU, Vissink A. **Periodontal inflamed surface area: quantifying inflammatory burden.** J Clin Periodontol. 2008;35(8):668-73. <https://pubmed.ncbi.nlm.nih.gov/18564145/>
74. Tonetti MS, D'Aiuto F, Nibali L, et al. **Treatment of periodontitis and endothelial function.** N Engl J Med. 2007;356(9):911-20. <https://pubmed.ncbi.nlm.nih.gov/17329698/>
75. D'Aiuto F, Parkar M, Andreou G, et al. **Periodontitis and systemic inflammation: control of the local infection is associated with a reduction in serum inflammatory markers.** J Dent Res. 2004;83(2):156-60. <https://pubmed.ncbi.nlm.nih.gov/14742655/>
76. Simpson TC, Clarkson JE, Worthington HV, et al. **Treatment of periodontitis for glycaemic control in people with diabetes mellitus.** Cochrane Database Syst Rev. 2022;4(4):CD004714. <https://pubmed.ncbi.nlm.nih.gov/35420698/>
77. Engebretson SP, Hyman LG, Michalowicz BS, et al. **The effect of nonsurgical periodontal therapy on hemoglobin A1c in persons with type 2 diabetes and chronic periodontitis: a randomized clinical trial (DPTT).** JAMA. 2013;310(23):2523-32. <https://pubmed.ncbi.nlm.nih.gov/24346989/>
78. Sanz M, Marco del Castillo A, Jepsen S, et al. **Periodontitis and cardiovascular diseases: consensus report.** J Clin Periodontol. 2020;47(3):268-88. <https://pubmed.ncbi.nlm.nih.gov/32011025/>
79. Michalowicz BS, Hodges JS, DiAngelis AJ, et al. **Treatment of periodontal disease and the risk of preterm birth.** N Engl J Med. 2006;355(18):1885-94. <https://pubmed.ncbi.nlm.nih.gov/17079762/>
80. Dominy SS, Lynch C, Ermini F, et al. **Porphyromonas gingivalis in Alzheimer's disease brains: evidence for disease causation and treatment with small-molecule inhibitors.** Sci Adv. 2019;5(1):eaau3333. <https://pubmed.ncbi.nlm.nih.gov/30746447/>

## 15. Classification · epidemiology · maintenance therapy

The basis for the `Staging` output, the relapse-prevention structure of maintenance therapy (SPT), and the `H0_TOOTH` hazard function.

81. Tonetti MS, Greenwell H, Kornman KS. **Staging and grading of periodontitis: framework and proposal of a new classification and case definition.** J Clin Periodontol. 2018;45 Suppl 20:S149-61. <https://pubmed.ncbi.nlm.nih.gov/29926495/>
82. Papapanou PN, Sanz M, Buduneli N, et al. **Periodontitis: consensus report of workgroup 2 of the 2017 World Workshop.** J Clin Periodontol. 2018;45 Suppl 20:S162-70. <https://pubmed.ncbi.nlm.nih.gov/29926490/>
83. Kassebaum NJ, Bernabé E, Dahiya M, Bhandari B, Murray CJ, Marcenes W. **Global burden of severe periodontitis in 1990-2010: a systematic review and meta-regression.** J Dent Res. 2014;93(11):1045-53. <https://pubmed.ncbi.nlm.nih.gov/25261053/>
84. Axelsson P, Nyström B, Lindhe J. **The long-term effect of a plaque control program on tooth mortality, caries and periodontal disease in adults. Results after 30 years of maintenance.** J Clin Periodontol. 2004;31(9):749-57. <https://pubmed.ncbi.nlm.nih.gov/15312097/>
85. Chambrone L, Chambrone D, Lima LA, Chambrone LA. **Predictors of tooth loss during long-term periodontal maintenance: a systematic review of observational studies.** J Clin Periodontol. 2010;37(7):675-84. <https://pubmed.ncbi.nlm.nih.gov/20528960/>
86. Lang NP, Adler R, Joss A, Nyman S. **Absence of bleeding on probing. An indicator of periodontal stability.** J Clin Periodontol. 1990;17(10):714-21. <https://pubmed.ncbi.nlm.nih.gov/2262585/>
87. Socransky SS, Haffajee AD, Goodson JM, Lindhe J. **New concepts of destructive periodontal disease.** J Clin Periodontol. 1984;11(1):21-32. <https://pubmed.ncbi.nlm.nih.gov/6584293/>
    (burst theory — the literature that corresponds directly to this model's prediction of site-level bistability)

## 16. Tissue regeneration

The basis for the `REGD` compartment and for `REGEN_FILL` (regeneration being the only route that reverses bone loss).

88. Nyman S, Lindhe J, Karring T, Rylander H. **New attachment following surgical treatment of human periodontal disease.** J Clin Periodontol. 1982;9(4):290-6. <https://pubmed.ncbi.nlm.nih.gov/6964676/>
89. Heijl L, Heden G, Svärdström G, Ostgren A. **Enamel matrix derivative (EMDOGAIN) in the treatment of intrabony periodontal defects.** J Clin Periodontol. 1997;24(9 Pt 2):705-14. <https://pubmed.ncbi.nlm.nih.gov/9310876/>
90. Nevins M, Giannobile WV, McGuire MK, et al. **Platelet-derived growth factor stimulates bone fill and rate of attachment level gain: results of a large multicenter randomized controlled trial.** J Periodontol. 2005;76(12):2205-15. <https://pubmed.ncbi.nlm.nih.gov/16332231/>

## 17. QSP methodology

91. Ribba B, Grimm HP, Agoram B, et al. **Methodologies for quantitative systems pharmacology (QSP) models: design and estimation.** CPT Pharmacometrics Syst Pharmacol. 2017;6(8):496-8. <https://pubmed.ncbi.nlm.nih.gov/28571119/>
92. Baron KT, Gastonguay MR. **Simulation from ODE-based population PK/PD and systems pharmacology models in R with mrgsolve.** J Pharmacokinet Pharmacodyn. 2015;42:S84-5. <https://mrgsolve.org/>
93. Cheng Y, Thalhauser CJ, Smithline S, et al. **QSP toolbox: computational implementation of integrated workflow components for deploying multi-scale mechanistic models.** AAPS J. 2017;19(4):1002-16. <https://pubmed.ncbi.nlm.nih.gov/28540623/>

---

## Quantitative calibration targets used directly in calibrating the model

The values below are the targets `PDT_calibration_report()` prints, and **not one of them was used as an
objective function.** The parameters were set from the mechanism and from independent measurements in the literature above, and
the clinical results were checked against them afterwards.

| Item | Literature value | Source |
|------|--------|------|
| SRP, change in PPD at 1-3 mm sites | -0.03 mm | Cobb 1996/2002 (46, 47) |
| SRP, change in CAL at 1-3 mm sites | -0.34 mm (a loss) | Cobb (46, 47) |
| SRP, reduction in PPD at 4-6 mm sites | 1.29 mm | Cobb (46, 47) |
| SRP, CAL gain at 4-6 mm sites | 0.55 mm | Cobb (46, 47) |
| SRP, reduction in PPD at ≥7 mm sites | 2.16 mm | Cobb (46, 47) |
| SRP, CAL gain at ≥7 mm sites | 1.19 mm | Cobb (46, 47) |
| Critical probing depth, SRP | 2.9 mm | Lindhe 1982 (48) |
| Critical probing depth, surgery | 4.2 mm | Lindhe 1982 (48) |
| Calculus removal rate at 3.0 / 4.5 / 6.5 mm | 83 / 61 / 32 % | Waerhaug 1978 (43) |
| Additional CAL gain with amoxicillin + metronidazole (deep sites) | 0.40-0.50 mm | Feres 2012 (53) |
| Additional CAL gain with SDD doxycycline | 0.30-0.40 mm | Caton 2000 (59) |
| SDD, reduction in GCF collagenase | 60-70 % | Golub 1998 (36) |
| Additional PPD reduction with local minocycline | 0.25-0.30 mm | Williams 2001 (57) |
| Additional PPD reduction with a chlorhexidine rinse | ~0.00 mm | Herrera 2020 (56) |
| AMY-101 alone, reduction in gingival inflammation at 28 days | 46 % | Hasturk 2021 (62) |
| Untreated alveolar bone loss | 0.15-0.40 mm/year | Löe 1986 (72) |
| Fall in HbA1c after periodontal treatment | 0.43 % | Simpson 2022 (76) |
| Fall in CRP after periodontal treatment | 0.4-0.9 mg/L | D'Aiuto 2004 (75) |
| Improvement in FMD at 6 months | +2.0 % | Tonetti 2007 (74) |
| FMD 24 hours after intensive treatment | transiently worse | Tonetti 2007 (74) |
| PISA, generalised stage III periodontitis | 1000-2000 mm² | Nesse 2008 (73) |
| Reduced PPD reduction in smokers | -0.40 to -0.50 mm | Labriola 2005 (69) |
| Distance from the plaque front to the bone crest | 1.0-2.0 mm | Waerhaug 1979 (33) |

---

## Disclaimer

This model is a qualitative to semi-quantitative QSP model for educational and hypothesis-generating purposes.
It rests on the published literature but has not been independently verified, and it must not be used directly in
clinical decision-making or in regulatory submission.
