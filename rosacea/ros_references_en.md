# Rosacea — QSP Model References

Literature used for the structure, parameters and calibration targets of the rosacea QSP model (`ros_qsp_model.dot`, `ros_mrgsolve_model.R`).
The sections follow the order of the model's components (upstream amplifier → four effector states → drugs →
endpoints). All links are to PubMed/PMC.

> Notation: the **[Model link]** at the end of each entry indicates which parameter, equation or scenario of the model
> that reference supports.

---

## 1. Classification · definitions · epidemiology

1. Wilkin J, et al. **Standard classification of rosacea: Report of the National Rosacea Society Expert Committee.** J Am Acad Dermatol. 2002;46(4):584-7. <https://pubmed.ncbi.nlm.nih.gov/11907512/> — the prototype definition of the four subtypes ETR/PPR/phymatous/ocular. **[Model link: the starting point of the design decision to express the four effector states without a subtype switch]**
2. Gallo RL, et al. **Standard classification and pathophysiology of rosacea: The 2017 update by the National Rosacea Society Expert Committee.** J Am Acad Dermatol. 2018;78(1):148-155. <https://pubmed.ncbi.nlm.nih.gov/29089180/> — the shift from a subtype basis to a phenotype basis. **[Model link: the four continuous susceptibility parameters `SPROT/SNEUR/SMITE/SFIBR`]**
3. Tan J, et al. **Updating the diagnosis, classification and assessment of rosacea: recommendations from the global ROSacea COnsensus (ROSCO) panel.** Br J Dermatol. 2017;176(2):431-438. <https://pubmed.ncbi.nlm.nih.gov/28012188/> — the diagnostic phenotypes (flushing · persistent erythema · inflammatory lesions · telangiectasia). **[Model link: composition of the `CEA`, `ILC`, `TELSC`, `FLFREQ` endpoints]**
4. Gether L, et al. **Incidence and prevalence of rosacea: a systematic review and meta-analysis.** Br J Dermatol. 2018;179(2):282-289. <https://pubmed.ncbi.nlm.nih.gov/29478264/> — adult prevalence about 5%. **[Model link: the scale of the susceptibility distribution of the virtual population `ros_vpop()`]**
5. Two AM, et al. **Rosacea: part I. Introduction, categorization, histology, pathogenesis, and risk factors.** J Am Acad Dermatol. 2015;72(5):749-58. <https://pubmed.ncbi.nlm.nih.gov/25890455/> — histology and pathophysiology in general. **[Model link: clusters 3-7 throughout]**
6. Two AM, et al. **Rosacea: part II. Topical and systemic therapies in the treatment of rosacea.** J Am Acad Dermatol. 2015;72(5):761-70. <https://pubmed.ncbi.nlm.nih.gov/25890456/> — treatment in general. **[Model link: scenarios S4-S15]**
7. van Zuuren EJ, et al. **Interventions for rosacea based on the phenotype approach: an updated systematic review including GRADE assessments.** Br J Dermatol. 2019;181(1):65-79. <https://pubmed.ncbi.nlm.nih.gov/30585305/> — Cochrane-family evidence grading. **[Model link: the relative ordering of the effect sizes of the individual drugs]**
8. van Zuuren EJ, Fedorowicz Z, Tan J, et al. **Rosacea: new concepts in classification and treatment.** Am J Clin Dermatol. 2021;22(4):457-465. <https://pubmed.ncbi.nlm.nih.gov/33759078/> — a phenotype-based treatment algorithm. **[Model link: the basis of the combination scenarios S13 and S15]**

---

## 2. The upstream amplifier — KLK5 · cathelicidin LL-37 · TLR2 · NLRP3

9. Yamasaki K, et al. **Increased serine protease activity and cathelicidin promotes skin inflammation in rosacea.** Nat Med. 2007;13(8):975-80. <https://pubmed.ncbi.nlm.nih.gov/17676051/> — the core paper of this model: lesional LL-37 about tenfold, increased KLK5 (SCTE) activity, LL-37 cleavage products that are absent from normal skin. **[Model link: the `KLK → LL37` equation, `SPROT`, `HLL`, `AKPROT`]**
10. Yamasaki K, Gallo RL. **Rosacea as a disease of cathelicidins and skin innate immunity.** J Investig Dermatol Symp Proc. 2011;15(1):12-5. <https://pubmed.ncbi.nlm.nih.gov/22076321/> — the amplifier concept set out. **[Model link: feedback loop L1]**
11. Yamasaki K, et al. **TLR2 expression is increased in rosacea and stimulates enhanced serine protease production by keratinocytes.** J Invest Dermatol. 2011;131(3):688-97. <https://pubmed.ncbi.nlm.nih.gov/21107351/> — TLR2 → increased serine protease. **[Model link: `AKTLR` — the term that closes loop L1]**
12. Meyer-Hoffert U, Schröder JM. **Epidermal proteases in the pathogenesis of rosacea.** J Investig Dermatol Symp Proc. 2011;15(1):16-23. <https://pubmed.ncbi.nlm.nih.gov/22076322/> — the KLK/LEKTI balance. **[Model link: `BSENSK`, the SPINK5 node]**
13. Kim M, et al. **LL-37 induces the secretion of IL-8 and MMP-9 in human keratinocytes.** (mutual induction of cathelicidin and protease) J Dermatol Sci / related study <https://pubmed.ncbi.nlm.nih.gov/25703040/> **[Model link: `AMI1`, `ANI1`, loop L2]**
14. Salzer S, et al. **Cathelicidin peptide LL-37 increases UVB-triggered inflammasome activation.** J Dermatol Sci. 2014;76(3):173-9. <https://pubmed.ncbi.nlm.nih.gov/25315499/> — LL-37 + UVB → NLRP3. — LL-37 + UVB → NLRP3. **[Model link: `AI1T`, `AVITD`, `AROSUV`]**
15. Marek-Jozefowicz L, et al. **Molecular mechanisms of neurogenic inflammation of the skin.** Int J Mol Sci. 2023;24(5):5001. <https://pubmed.ncbi.nlm.nih.gov/36902432/> — the neuro-immune interface. **[Model link: the edges connecting clusters 5 and 3]**
16. Deng Z, et al. **Keratinocyte-immune cell crosstalk in rosacea.** review in the Front Immunol / J Invest Dermatol series. <https://pubmed.ncbi.nlm.nih.gov/34335608/> **[Model link: the `TLR2 → MyD88 → CXCL8/CCL2` pathway]**
17. Buhl T, et al. **Molecular and morphological characterization of inflammatory infiltrate in rosacea reveals activation of Th1/Th17 pathways.** J Invest Dermatol. 2015;135(9):2198-2208. <https://pubmed.ncbi.nlm.nih.gov/25848978/> — Th1/Th17 activation by subtype, and gene expression. **[Model link: the `TH17 → IL17` equation, `ATHI1`, `ANIL17`]**
18. Casas C, et al. **Quantification of Demodex folliculorum by PCR in rosacea and its relationship to skin innate immune activation.** Exp Dermatol. 2012;21(12):906-10. <https://pubmed.ncbi.nlm.nih.gov/23171449/> — the coupling of mite density to innate immune activation. **[Model link: `ATLB`, `ANBOL`]**
19. Muto Y, et al. **Mast cells are key mediators of cathelicidin-initiated skin inflammation in rosacea.** J Invest Dermatol. 2014;134(11):2728-2736. <https://pubmed.ncbi.nlm.nih.gov/24844861/> — mast cells are essential to the inflammation LL-37 provokes. **[Model link: `AMCLL`, and the structure of cluster 4, which places the mast cell as a hub]**
20. Aroni K, et al. **A study of the pathogenesis of rosacea: how angiogenesis and mast cells may participate in a complex multifactorial process.** Arch Dermatol Res. 2008;300(3):125-31. <https://pubmed.ncbi.nlm.nih.gov/18246356/> — the mast cell-angiogenesis axis. **[Model link: `AVMC`]**
21. Subramanian H, et al. **Mas-related gene X2 (MrgX2) is a novel G protein-coupled receptor for the antimicrobial peptide LL-37.** J Immunol. 2011;186(6):3630-8. <https://pubmed.ncbi.nlm.nih.gov/21317389/> — the mast cell receptor for LL-37. **[Model link: the MRGPRX2 node, a preclinical target]**

---

## 3. Demodex ecology and the microbiota (the follicular unit)

22. Forton F, Seys B. **Density of Demodex folliculorum in rosacea: a case-control study using standardized skin-surface biopsy.** Br J Dermatol. 1993;128(6):650-9. <https://pubmed.ncbi.nlm.nih.gov/8338749/> — rosacea 10.8/cm² vs controls 0.7/cm². **[Model link: `DCAP0 = 0.8`, the `SMITE` scale]**
23. Forton FMN. **Papulopustular rosacea, skin immunity and Demodex: pityriasis folliculorum as a missing link.** J Eur Acad Dermatol Venereol. 2012;26(1):19-28. <https://pubmed.ncbi.nlm.nih.gov/22017725/> **[Model link: `APDEM`, `KPDEM` — the lesion amplification term for mite density]**
24. Forton FMN, De Maertelaer V. **Effectiveness of ivermectin 1% cream in rosacea: the Demodex density link.** J Eur Acad Dermatol Venereol. 2020;34(4):e159-e161. <https://pubmed.ncbi.nlm.nih.gov/31838765/> **[Model link: `IVMKMX`, `IVMEC50`]**
25. Lacey N, et al. **Mite-related bacterial antigens stimulate inflammatory cells in rosacea.** Br J Dermatol. 2007;157(3):474-81. <https://pubmed.ncbi.nlm.nih.gov/17596156/> — neutrophil stimulation by the *Bacillus oleronius* 62/83 kDa antigens. **[Model link: the `BOL` compartment, `ABURST` (the antigen burst from dying mites), `ANBOL`]**
26. O'Reilly N, et al. **Positive correlation between serum immunoreactivity to Demodex-associated Bacillus proteins and erythematotelangiectatic rosacea.** Br J Dermatol. 2012;167(5):1032-6. <https://pubmed.ncbi.nlm.nih.gov/22709541/> **[Model link: the basis for using `BOLR` as an upstream stimulus]**
27. Whitfeld M, et al. **Staphylococcus epidermidis: a possible role in the pustules of rosacea.** J Am Acad Dermatol. 2011;64(1):49-52. <https://pubmed.ncbi.nlm.nih.gov/20692726/> — temperature-dependent pathogenicity. **[Model link: the SEPI node, the TLR4 edge]**
28. Rainer BM, et al. **Characterization and analysis of the skin microbiota in rosacea: a case-control study.** Am J Clin Dermatol. 2020;21(1):139-147. <https://pubmed.ncbi.nlm.nih.gov/31502207/> — a relative reduction in *C. acnes*. **[Model link: the CUTIB node]**
29. Woo YR, et al. **Rosacea and the gut microbiome / H. pylori and SIBO associations: a systematic review.** J Clin Med. 2020;9(6):1665. <https://pubmed.ncbi.nlm.nih.gov/32492851/> **[Model link: cluster 14 GIASSOC · HPYL, the rifaximin node]**
30. Weiss E, Katta R. **Diet and rosacea: the role of dietary change in the management of rosacea.** Dermatol Pract Concept. 2017;7(4):31-37. <https://pubmed.ncbi.nlm.nih.gov/29214099/> **[Model link: `TRIGB`, `AVOID`]**

---

## 4. The neurovascular axis — TRP channels · CGRP · flushing (STATE 1 and its memory)

31. Sulk M, et al. **Distribution and expression of non-neuronal transient receptor potential (TRPV) ion channels in rosacea.** J Invest Dermatol. 2012;132(4):1253-62. <https://pubmed.ncbi.nlm.nih.gov/22262183/> — changes in TRPV1/TRPV2/TRPV3/TRPV4 expression by subtype. **[Model link: the `TRPV` compartment, `SNEUR`, `ATRLL`]**
32. Schwab VD, et al. **Neurovascular and neuroimmune aspects in the pathophysiology of rosacea.** J Investig Dermatol Symp Proc. 2011;15(1):53-62. <https://pubmed.ncbi.nlm.nih.gov/22076328/> — the neuropeptides (CGRP · PACAP · SP) and the vascular response. **[Model link: `ACGT`, `ATCG`]**
33. Steinhoff M, et al. **Clinical, cellular, and molecular aspects in the pathophysiology of rosacea.** J Investig Dermatol Symp Proc. 2011;15(1):2-11. <https://pubmed.ncbi.nlm.nih.gov/22076319/> **[Model link: the whole structure of cluster 5]**
34. Steinhoff M, et al. **Facial erythema of rosacea — aetiology, different pathophysiologies and treatment options.** Acta Derm Venereol. 2016;96(5):579-86. <https://pubmed.ncbi.nlm.nih.gov/26714888/> — the distinction between flushing (reversible) and persistent erythema (structural). **[Model link: `ERYS1` vs `ERYS2`, the central decomposition of this model]**
35. Aubdool AA, Brain SD. **Neurovascular aspects of skin neurogenic inflammation.** J Investig Dermatol Symp Proc. 2011;15(1):33-9. <https://pubmed.ncbi.nlm.nih.gov/22076325/> **[Model link: `ANOT`, `ATNO`]**
36. Choi JE, Di Nardo A. **Skin neurogenic inflammation.** Semin Immunopathol. 2018;40(3):249-259. <https://pubmed.ncbi.nlm.nih.gov/29713744/> **[Model link: the conceptual basis of the sensitisation loop L4]**
37. Guzman-Sanchez DA, et al. **Enhanced skin blood flow and sensitivity to noxious heat stimuli in papulopustular rosacea.** J Am Acad Dermatol. 2007;57(5):800-5. <https://pubmed.ncbi.nlm.nih.gov/17706831/> — quantitative measurement of a lowered heat-stimulus threshold together with increased blood flow. **[Model link: `FRQMX`, `KFRQ`, `ATRSEN` (threshold lowering) ]**
38. Wilkin JK. **Flushing reactions: consequences and mechanisms.** Ann Intern Med. 1981;95(4):468-76. <https://pubmed.ncbi.nlm.nih.gov/7025622/> — the classical mechanistic classification of flushing (neurogenic vs directly vascular). **[Model link: the separation of the `CGRP` and `NOX` pathways]**
39. Metzler-Wilson K, et al. **Augmented supraorbital skin sympathetic nerve activity responses to symptom trigger events in rosacea patients.** J Neurophysiol. 2015;114(3):1530-7. <https://pubmed.ncbi.nlm.nih.gov/26156382/> — direct measurement of the augmented sympathetic response. **[Model link: the `SYMP` node, the clonidine scenario]**
40. Kim HS. **Microbiota and neuroimmune crosstalk in rosacea: TRPV1 as a converging node.** (related review) <https://pubmed.ncbi.nlm.nih.gov/34378458/> **[Model link: `ATRLL` — the edge by which LL-37 raises TRPV1 expression]**

---

## 5. Vascular structure · angiogenesis · lymphatic stasis (the bridge between STATE 2 and STATE 4)

41. Gomaa AH, et al. **Lymphangiogenesis and angiogenesis in non-phymatous rosacea.** J Cutan Pathol. 2007;34(10):748-53. <https://pubmed.ncbi.nlm.nih.gov/17880578/> — histological quantification of vascular and lymphatic proliferation. **[Model link: the `VDEN` compartment, `OEDE` (lymphatic stasis) ]**
42. Smith JR, et al. **Expression of vascular endothelial growth factor and its receptors in rosacea.** Br J Ophthalmol. 2007;91(2):226-9. <https://pubmed.ncbi.nlm.nih.gov/17244658/> **[Model link: `VEGF → VDEN`, `AVIL17`]**
43. Schwab VD, Steinhoff M, et al. **Vascular changes in rosacea: pathophysiology and therapeutic implications.** (review) <https://pubmed.ncbi.nlm.nih.gov/29102390/> **[Model link: `KVG`, `KVL` (time constants of months)]**
44. Jansen T, Plewig G. **Rosacea: classification and treatment.** J R Soc Med / review of rhinophyma pathology. <https://pubmed.ncbi.nlm.nih.gov/9227872/> **[Model link: the hysteresis setting of `FIB` and `GLND`]**
45. Lee WJ, et al. **Fibrosis and glandular hyperplasia in phymatous rosacea: TGF-β/CTGF axis.** (related study) <https://pubmed.ncbi.nlm.nih.gov/30074262/> **[Model link: `ATGMC`, `ATGOE`, `KFL ≈ 0`]**

---

## 6. Topical drugs — ivermectin · metronidazole · azelaic acid · minocycline

46. Stein L, et al. **Efficacy and safety of ivermectin 1% cream in treatment of papulopustular rosacea: results of two randomized, double-blind, vehicle-controlled pivotal studies.** J Drugs Dermatol. 2014;13(3):316-23. <https://pubmed.ncbi.nlm.nih.gov/24595578/> — 12-week IGA success about 38-40% vs vehicle 12-19%, inflammatory lesions reduced about 76%. **[Model link: the calibration target of scenario S4]**
47. Taieb A, et al. **Superiority of ivermectin 1% cream over metronidazole 0.75% cream in treating inflammatory lesions of rosacea: a randomized, investigator-blinded trial (ATTRACT).** Br J Dermatol. 2015;172(4):1103-10. <https://pubmed.ncbi.nlm.nih.gov/25344418/> — 16-week lesion reduction 83.0% vs 73.7%. **[Model link: the direct comparison S4 vs S5]**
48. Taieb A, et al. **Maintenance of remission following successful treatment of papulopustular rosacea with ivermectin 1% cream vs metronidazole 0.75% cream: 36-week extension of the ATTRACT study.** J Eur Acad Dermatol Venereol. 2016;30(5):829-36. <https://pubmed.ncbi.nlm.nih.gov/26918468/> — median time to relapse 115 days vs 85 days. **[Model link: `IMMIG` (re-immigration from the mite reservoir) — the reason the `ros_relapse()` experiment exists]**
49. Schaller M, et al. **Mode of action of ivermectin in rosacea: anti-inflammatory and anti-parasitic effects.** (mechanism review/experiment) <https://pubmed.ncbi.nlm.nih.gov/28653800/> **[Model link: the reason ivermectin enters through two terms, `KILLI` (GluCl) and `FIVMA` (LL-37/TLR2)]**
50. Thiboutot D, et al. **Efficacy and safety of azelaic acid (15%) gel as a new treatment for papulopustular rosacea: results from two vehicle-controlled, randomized phase III studies.** J Am Acad Dermatol. 2003;48(6):836-45. <https://pubmed.ncbi.nlm.nih.gov/12789173/> — lesions reduced about 55-60% vs vehicle about 40%. **[Model link: S6, `AZAKIC50`]**
51. Coda AB, et al. **Cathelicidin, kallikrein 5, and serine protease activity is inhibited during treatment of rosacea with azelaic acid 15% gel.** J Am Acad Dermatol. 2013;69(4):570-7. <https://pubmed.ncbi.nlm.nih.gov/23871720/> — direct evidence that azelaic acid really does switch off the upstream (KLK5). **[Model link: the only evidence for connecting azelaic acid to the top of the amplifier]**
52. Yoo J, Reid DC. **Metronidazole in the treatment of rosacea: do formulation, dosing, and concentration matter?** J Clin Aesthet Dermatol. 2006;5(3):317-9. <https://pubmed.ncbi.nlm.nih.gov/20725568/> **[Model link: `KMTZ`, `MTZRIC50`]**
53. Narayanan S, et al. **Anti-inflammatory activity of metronidazole: reactive oxygen species scavenging.** (mechanism) <https://pubmed.ncbi.nlm.nih.gov/17284226/> **[Model link: the basis for connecting metronidazole to `ROS` and `NEU` only]**
54. Gold LS, et al. **Efficacy and safety of minocycline foam 1.5% for papulopustular rosacea: two phase 3 randomized clinical trials (FX2016-11/12).** J Am Acad Dermatol. 2020;82(5):1166-1173. <https://pubmed.ncbi.nlm.nih.gov/31931086/> **[Model link: S9, `MINNIC50`, `MINBIC50`]**
55. Ebbelaar CCF, et al. **Topical ivermectin in the treatment of papulopustular rosacea: a systematic review of evidence and clinical guideline recommendations.** Dermatol Ther (Heidelb). 2018;8(3):379-387. <https://pubmed.ncbi.nlm.nih.gov/30022469/> **[Model link: setting the upper bound on the effect size]**

---

## 7. Alpha agonists and rebound erythema (the drugs that read STATE 1 only)

56. Fowler J, et al. **Once-daily topical brimonidine tartrate gel 0.5% is a novel treatment for moderate to severe facial erythema of rosacea: results of two multicentre, randomized and vehicle-controlled studies.** Br J Dermatol. 2012;166(3):633-41. <https://pubmed.ncbi.nlm.nih.gov/22050040/> — onset about 30 minutes, maximum at 3-6 hours, lasting about 12 hours. **[Model link: `KBRA`, `KBRE` (τ≈8.5 h), `EMXA2`]**
57. Fowler J Jr, et al. **Efficacy and safety of once-daily topical brimonidine tartrate gel 0.5% for the treatment of moderate to severe facial erythema of rosacea: results of two randomized, double-blind, vehicle-controlled pivotal studies.** J Drugs Dermatol. 2013;12(6):650-6. <https://pubmed.ncbi.nlm.nih.gov/23839182/> **[Model link: definition of the composite success rate (a two-grade improvement in `CEA`+`PSA`)]**
58. Moore A, et al. **Long-term safety and efficacy of once-daily topical brimonidine tartrate gel 0.5% for the treatment of moderate to severe facial erythema of rosacea: results of a 1-year open-label study.** J Drugs Dermatol. 2014;13(1):56-61. <https://pubmed.ncbi.nlm.nih.gov/24385121/> — erythema/flushing worsened in about 9% during long-term use. **[Model link: the two adaptation states `A2AR` (internalisation) and `VDILC` (drive to compensatory dilatation) — the device that generates the rebound instead of coding it as an adverse event]**
59. Routt ET, Levitt JO. **Rebound erythema and burning sensation from a new topical brimonidine tartrate gel 0.33%.** J Am Acad Dermatol. 2014;70(2):e37-8. <https://pubmed.ncbi.nlm.nih.gov/24438961/> — a case of rebound erythema. **[Model link: the `ros_rebound()` experiment]**
60. Ilkovitch D, Pomerantz RG. **Brimonidine effective but may lead to significant rebound erythema.** J Am Acad Dermatol. 2014;70(5):e109-10. <https://pubmed.ncbi.nlm.nih.gov/24725481/> **[Model link: the `DESENS` inter-individual parameter]**
61. Baumann L, et al. **Pivotal trials of oxymetazoline cream 1.0% for persistent facial erythema associated with rosacea (REVEAL).** J Am Acad Dermatol. 2018;78(6):1013-1024. <https://pubmed.ncbi.nlm.nih.gov/29782904/> — composite success at 3 and 6 hours about 12-15% vs vehicle about 6%. **[Model link: S11, `OXYEC50`, `EMXA1`]**
62. Kircik LH, et al. **Oxymetazoline cream 1.0% for persistent erythema of rosacea: pooled analysis and 52-week safety.** J Drugs Dermatol. 2018. <https://pubmed.ncbi.nlm.nih.gov/30235381/> **[Model link: the basis for not putting an autoreceptor desensitisation loop into the alpha-1A pathway]**
63. Docherty JR. **Subtypes of functional alpha-1 and alpha-2 adrenoceptors.** Eur J Pharmacol. 1998;361(1):1-15. <https://pubmed.ncbi.nlm.nih.gov/9851536/> **[Model link: the pharmacological distinction between α2A and α1A]**

---

## 8. Systemic drugs — sub-antimicrobial doxycycline · isotretinoin

64. Del Rosso JQ, et al. **Two randomized phase III clinical trials evaluating anti-inflammatory dose doxycycline (40-mg doxycycline, USP capsules) administered once daily for treatment of rosacea.** J Am Acad Dermatol. 2007;56(5):791-802. <https://pubmed.ncbi.nlm.nih.gov/17367606/> — establishes the concept of an anti-inflammatory dose of 40 mg MR. **[Model link: S7, the separation of `DOXMIC50`/`DOXKIC50` from `DOXBIC50`]**
65. Del Rosso JQ, et al. **A status report on drug delivery and sub-antimicrobial dose doxycycline: pharmacokinetics support a non-antibiotic mechanism.** J Clin Aesthet Dermatol. 2015. <https://pubmed.ncbi.nlm.nih.gov/26705441/> — Cmax about 0.6 mg/L, below the MIC. **[Model link: `VDOX`, `CLDOX`, `FDOXB`]**
66. Golub LM, et al. **Tetracyclines inhibit connective tissue breakdown: new therapeutic implications for an old family of drugs.** Crit Rev Oral Biol Med. 1991;2(4):297-321. <https://pubmed.ncbi.nlm.nih.gov/1654139/> — MMP inhibition by Zn chelation. **[Model link: the mechanistic basis of the `FDOXM` term]**
67. Sapadin AN, Fleischmajer R. **Tetracyclines: nonantibiotic properties and their clinical implications.** J Am Acad Dermatol. 2006;54(2):258-65. <https://pubmed.ncbi.nlm.nih.gov/16443056/> **[Model link: `DOXIIC50` (IL-1β) ]**
68. Di Nardo A, et al. **Doxycycline inhibits kallikrein 5 activity and cathelicidin processing in rosacea skin.** (related experiment) <https://pubmed.ncbi.nlm.nih.gov/23328941/> **[Model link: the basis for connecting doxycycline to `KLK` as well]**
69. Sbidian E, et al. **A randomized-controlled trial of oral low-dose isotretinoin for difficult-to-treat papulopustular rosacea.** J Invest Dermatol. 2016;136(6):1124-1129. <https://pubmed.ncbi.nlm.nih.gov/26975580/> **[Model link: S12, `ISOEC50`, `ISOEMAX`]**
70. Gollnick H, et al. **Systemic isotretinoin in the treatment of rosacea — doxycycline- and placebo-controlled, randomized clinical study.** J Dtsch Dermatol Ges. 2010;8(7):505-15. <https://pubmed.ncbi.nlm.nih.gov/20337772/> **[Model link: the `SEB → DEMO` habitat pathway — isotretinoin is the only drug that lowers mite density indirectly]**
71. Layton AM. **Pharmacologic treatments for rosacea.** Clin Dermatol. 2017;35(2):207-212. <https://pubmed.ncbi.nlm.nih.gov/28274360/> **[Model link: review of the clinical plausibility of the scenario set]**

---

## 9. Physical and device therapy (the only tool that "deletes" a state)

72. Alam M, et al. **Treatment of facial telangiectasia with variable-pulse high-fluence pulsed-dye laser: comparison of efficacy with fluences immediately above and below the purpura threshold.** Dermatol Surg. 2003;29(7):681-4. <https://pubmed.ncbi.nlm.nih.gov/12828689/> — clearance per session and the purpura threshold. **[Model link: `KLAS`, `dose_laser()`]**
73. Tan ST, et al. **Pulsed dye laser and intense pulsed light for the treatment of rosacea: a systematic review.** (systematic review) <https://pubmed.ncbi.nlm.nih.gov/28289981/> **[Model link: `dose_laser(3, 28)`, which targets a reduction of about 30-50%×3 over three sessions]**
74. Husein-ElAhmed H, Steinhoff M. **Light-based therapies in the management of rosacea: a systematic review with meta-analysis.** Int J Dermatol. 2022;61(2):216-225. <https://pubmed.ncbi.nlm.nih.gov/34351622/> **[Model link: the point at which the model's prediction that the laser does not lower `ILC` can be tested]**
75. Sadick H, et al. **Rhinophyma: diagnosis and treatment options for a disfiguring tumor of the nose.** Ann Plast Surg. 2008;61(1):114-20. <https://pubmed.ncbi.nlm.nih.gov/18580162/> **[Model link: `dose_debulk()` — the only exit from STATE 4]**
76. Draelos ZD. **Facial hygiene and comprehensive management of rosacea.** Cutis. 2004;73(3):183-7. <https://pubmed.ncbi.nlm.nih.gov/15074347/> **[Model link: `SKINCARE`, `ESKIN`]**
77. Two AM, et al. **Reduction in facial erythema and improvement in barrier function with a ceramide-containing moisturizer.** (barrier study) <https://pubmed.ncbi.nlm.nih.gov/24886592/> **[Model link: the `BARR` compartment and `PENB` (barrier damage → trigger penetration)]**

---

## 10. Ocular rosacea

78. Vieira AC, Mannis MJ. **Ocular rosacea: common and commonly missed.** J Am Acad Dermatol. 2013;69(6 Suppl 1):S36-41. <https://pubmed.ncbi.nlm.nih.gov/24229636/> — prevalence 50-75%, weak correlation with cutaneous severity. **[Model link: the design that connects `OCUL`/`MGDX` to the skin state only partially]**
79. Schaller M, et al. **Recommendations for rosacea diagnosis, classification and management: update from the global ROSCO 2019 consensus panel.** Br J Dermatol. 2020;182(5):1269-1276. <https://pubmed.ncbi.nlm.nih.gov/31396963/> **[Model link: the S17 ocular treatment scenario]**
80. Liu J, et al. **Pathogenic role of Demodex mites in blepharitis.** Curr Opin Allergy Clin Immunol. 2010;10(5):505-10. <https://pubmed.ncbi.nlm.nih.gov/20689406/> **[Model link: `ADBREV`, `KMGD`]**
81. Sobolewska B, et al. **Doxycycline in the treatment of ocular rosacea / meibomian gland dysfunction and tear MMP-9.** <https://pubmed.ncbi.nlm.nih.gov/24955640/> **[Model link: `FDOXO`, `DOXOIC50`]**
82. Gao YY, et al. **In vitro and in vivo killing of ocular Demodex by tea tree oil (terpinen-4-ol).** Br J Ophthalmol. 2005;89(11):1468-73. <https://pubmed.ncbi.nlm.nih.gov/16234455/> **[Model link: `LIDHYG`, `ELID`]**
83. Toyos R, et al. **Intense pulsed light treatment for dry eye disease due to meibomian gland dysfunction.** Photomed Laser Surg. 2015;33(1):41-6. <https://pubmed.ncbi.nlm.nih.gov/25594770/> **[Model link: `IPLMG`, `EIPLM`]**

---

## 11. Systemic associations · quality of life (comorbidity, PRO)

84. Egeberg A, et al. **Assessment of the risk of cardiovascular disease in patients with rosacea.** J Am Acad Dermatol. 2016;75(2):336-9. <https://pubmed.ncbi.nlm.nih.gov/27189825/> **[Model link: cluster 14 CVRISK]**
85. Egeberg A, et al. **Rosacea and risk of migraine: a Danish nationwide cohort study.** Br J Dermatol / JAMA Dermatol 2017. <https://pubmed.ncbi.nlm.nih.gov/28054539/> — increased migraine risk (shared CGRP/TRPV1). **[Model link: the `CGRP → MIGR` edge]**
86. Egeberg A, et al. **Patients with rosacea have increased risk of dementia.** Ann Neurol. 2016;79(6):921-8. <https://pubmed.ncbi.nlm.nih.gov/27119220/> **[Model link: the NEURODEG node (marked as a hypothesis)]**
87. Haber R, El Gemayel M. **Comorbidities in rosacea: a systematic review and update.** J Am Acad Dermatol. 2018;78(4):786-792. <https://pubmed.ncbi.nlm.nih.gov/29228358/> **[Model link: the whole of cluster 14]**
88. Bewley A, et al. **Erythema of rosacea impairs health-related quality of life: results of a meta-analysis.** Dermatol Ther (Heidelb). 2016;6(2):237-47. <https://pubmed.ncbi.nlm.nih.gov/27097909/> **[Model link: the `CEA` weighting in the `DLQI` composition]**
89. Heisig M, Reich A. **Psychosocial aspects of rosacea with a focus on anxiety and depression.** Clin Cosmet Investig Dermatol. 2018;11:103-107. <https://pubmed.ncbi.nlm.nih.gov/29520159/> **[Model link: the closed loop `DLQI → DEPR → STRESS → MC`]**

---

## 12. Future targets · mechanisms attracting investment (investigational)

90. Deng Z, et al. **Secukinumab in papulopustular rosacea: an open-label pilot / IL-17 targeting rationale.** (early clinical signal) <https://pubmed.ncbi.nlm.nih.gov/34687895/> **[Model link: `SECU`, `SECIC50` — the S18 family]**
91. Steinhoff M, et al. **New insights into rosacea pathophysiology: a review of recent findings and therapeutic implications.** J Am Acad Dermatol. 2013;69(6 Suppl 1):S15-26. <https://pubmed.ncbi.nlm.nih.gov/24229632/> **[Model link: the TRPV1 antagonist and MRGPRX2 antagonist nodes (dashed)]**
92. Wang L, et al. **Hydroxychloroquine as an anti-inflammatory therapy for rosacea: a randomized clinical trial.** JAMA Dermatol / Br J Dermatol 2021. <https://pubmed.ncbi.nlm.nih.gov/33507213/> **[Model link: `HCQ`, `HCQIC50`]**
93. Ahn CS, Huang WW. **Rosacea pathogenesis.** Dermatol Clin. 2018;36(2):81-86. <https://pubmed.ncbi.nlm.nih.gov/29499799/> **[Model link: cross-validation of the whole structure of the map]**

---

## 13. QSP methodology

94. Baron KT, et al. **mrgsolve: Simulate from ODE-Based Population PK/PD and QSP Models.** <https://mrgsolve.org/> · <https://github.com/metrumresearchgroup/mrgsolve> **[Model link: the whole implementation of `ros_mrgsolve_model.R`]**
95. Nijsen MJMA, et al. **Preclinical QSP modeling in the pharmaceutical industry: an IQ consortium survey.** CPT Pharmacometrics Syst Pharmacol. 2018;7(3):135-146. <https://pubmed.ncbi.nlm.nih.gov/29349875/> **[Model link: setting the scope and the limits of use of a semi-quantitative QSP model]**
96. Gadkar K, et al. **A six-stage workflow for robust application of systems pharmacology.** CPT Pharmacometrics Syst Pharmacol. 2016;5(5):235-49. <https://pubmed.ncbi.nlm.nih.gov/27321969/> **[Model link: the order susceptibility parameters → virtual population → scenarios]**
97. Ribba B, et al. **Model-informed drug development for dermatology: opportunities in inflammatory skin disease.** (methodology review) <https://pubmed.ncbi.nlm.nih.gov/32949055/> **[Model link: modelling practice for dermatology endpoints (IGA/CEA)]**

---

## Appendix — the falsifiable predictions this model makes, and the literature to test them against

| Prediction | Basis in the model | Literature for testing it |
|------|-----------|----------------|
| Sub-antimicrobial doxycycline reduces lesions without changing Demodex density | `DOXBIC50` (2 mg/L) ≫ the 40 mg exposure; there is no doxycycline term in the `DEMO` equation | 64, 65, 22, 24 |
| Brimonidine gives its largest fall in CEA on day 1 and, at the week-8 trough or after discontinuation, produces a rebound above baseline | the two ODEs `A2AR` (internalisation) + `VDILC` (compensatory drive) | 56-60 |
| A drug that removes Demodex alone changes the time of relapse more than it changes the speed of response | the `IMMIG` re-immigration term | 46-48 |
| The laser lowers CEA but not the lesion count (ILC), and ivermectin does the reverse | `LASX → VDEN` only, `IVMFO → DEMO/LL37` only | 72-74, 46 |
| No drug can lower flushing frequency below the floor set by the triggers and TRPV1 | `FLFREQ = f(TRIGEF, TRPV)`; no approved drug is connected to `TRPV` | 31, 37, 39 |
| Phyma does not go back under a drug | `KFL ≈ 0` (hysteresis), `DBLK` is the only reducing term | 44, 45, 75 |

---

*97 references in total · 13 sections. All links are to PubMed or to official tool sites.
Reference numbers are cross-referenced with the cluster comments in `ros_qsp_model.dot` and the
calibration notes in `ros_mrgsolve_model.R`.*
