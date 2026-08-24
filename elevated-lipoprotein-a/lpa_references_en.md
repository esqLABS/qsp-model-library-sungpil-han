# Elevated Lipoprotein(a) — References
# Elevated Lipoprotein(a) — Reference List

This list sets out, section by section, the basis for the structure and parameters of
`lpa_qsp_model.dot` / `lpa_mrgsolve_model_en.R`. The PMID link on each entry goes to PubMed.
For the reference each equation of the model rests on, see the map in §13.

---

## 1. Discovery · structure · evolution

1. Berg K. **A new serum type system in man — the Lp system.** *Acta Pathol Microbiol Scand.* 1963;59:369-82. — the first description of Lp(a). [PMID 14064818](https://pubmed.ncbi.nlm.nih.gov/14064818/)
2. McLean JW, Tomlinson JE, Kuang WJ, et al. **cDNA sequence of human apolipoprotein(a) is homologous to plasminogen.** *Nature.* 1987;330(6144):132-7. — the paper showing apo(a) to be a homologue of plasminogen. The basis of the whole of ARM 3 in the model. [PMID 3670400](https://pubmed.ncbi.nlm.nih.gov/3670400/)
3. Utermann G. **The mysteries of lipoprotein(a).** *Science.* 1989;246(4932):904-10. [PMID 2530631](https://pubmed.ncbi.nlm.nih.gov/2530631/)
4. Lawn RM, Schwartz K, Patthy L. **Convergent evolution of apolipoprotein(a) in primates and hedgehog.** *Proc Natl Acad Sci USA.* 1997;94(22):11992-7. — present only in Old World primates and the hedgehog. [PMID 9342350](https://pubmed.ncbi.nlm.nih.gov/9342350/)
5. Koschinsky ML, Marcovina SM. **Structure-function relationships in apolipoprotein(a): insights into lipoprotein(a) assembly and pathogenicity.** *Curr Opin Lipidol.* 2004;15(2):167-74. [PMID 15017359](https://pubmed.ncbi.nlm.nih.gov/15017359/)
6. Schmidt K, Noureen A, Kronenberg F, Utermann G. **Structure, function, and genetics of lipoprotein(a).** *J Lipid Res.* 2016;57(8):1339-59. — comprehensive review. [PMID 27074913](https://pubmed.ncbi.nlm.nih.gov/27074913/)

## 2. Genetics · the KIV-2 copy-number variant

7. Utermann G, Menzel HJ, Kraft HG, et al. **Lp(a) glycoprotein phenotypes. Inheritance and relation to Lp(a)-lipoprotein concentrations in plasma.** *J Clin Invest.* 1987;80(2):458-65. — the inverse correlation between isoform size and plasma concentration. The basis of the `SECEFF` term in the model. [PMID 2956279](https://pubmed.ncbi.nlm.nih.gov/2956279/)
8. Boerwinkle E, Leffert CC, Lin J, et al. **Apolipoprotein(a) gene accounts for greater than 90% of the variation in plasma lipoprotein(a) concentrations.** *J Clin Invest.* 1992;90(1):52-60. — the reason the model makes Lp(a) a production-determined variable. [PMID 1386087](https://pubmed.ncbi.nlm.nih.gov/1386087/)
9. Kraft HG, Köchl S, Menzel HJ, Sandholzer C, Utermann G. **The apolipoprotein(a) gene: a transcribed hypervariable locus controlling plasma lipoprotein(a) concentration.** *Hum Genet.* 1992;90(3):220-30. [PMID 1487235](https://pubmed.ncbi.nlm.nih.gov/1487235/)
10. Clarke R, Peden JF, Hopewell JC, et al. **Genetic variants associated with Lp(a) lipoprotein level and coronary disease.** *N Engl J Med.* 2009;361(26):2518-28. — rs10455872 · rs3798220. [PMID 20032323](https://pubmed.ncbi.nlm.nih.gov/20032323/)
11. Kamstrup PR, Tybjaerg-Hansen A, Steffensen R, Nordestgaard BG. **Genetically elevated lipoprotein(a) and increased risk of myocardial infarction.** *JAMA.* 2009;301(22):2331-9. — Mendelian randomisation. [PMID 19509380](https://pubmed.ncbi.nlm.nih.gov/19509380/)
12. Mukamel RE, Handsaker RE, Sherman MA, et al. **Protein-coding repeat polymorphisms strongly modulate plasma protein levels.** *Science.* 2021;373(6562):1499-1505. — the quantitative effect of KIV-2 copy number. [PMID 34554798](https://pubmed.ncbi.nlm.nih.gov/34554798/)
13. Coassin S, Kronenberg F. **Lipoprotein(a) beyond the kringle IV repeat polymorphism: from genotype to phenotype.** *Atherosclerosis.* 2022;349:17-35. [PMID 35606073](https://pubmed.ncbi.nlm.nih.gov/35606073/)
14. Trinder M, Uddin MM, Finneran P, Aragam KG, Natarajan P. **Clinical utility of lipoprotein(a) and LPA genetic risk score in risk prediction of incident atherosclerotic cardiovascular disease.** *JAMA Cardiol.* 2021;6(3):287-95. [PMID 33021622](https://pubmed.ncbi.nlm.nih.gov/33021622/)
15. Guan W, Cao J, Steffen BT, et al. **Race is a key variable in assigning lipoprotein(a) cutoff values for coronary heart disease risk assessment: the MESA study.** *Arterioscler Thromb Vasc Biol.* 2015;35(4):996-1001. — the model's `FANC` (race) covariate. [PMID 25810300](https://pubmed.ncbi.nlm.nih.gov/25810300/)
16. Patel AP, Wang M, Pirruccello JP, et al. **Lp(a) (lipoprotein[a]) concentrations and incident atherosclerotic cardiovascular disease: new insights from a large national biobank.** *Arterioscler Thromb Vasc Biol.* 2021;41(1):465-74. [PMID 33115266](https://pubmed.ncbi.nlm.nih.gov/33115266/)

## 3. Biosynthesis · secretion efficiency · assembly

17. White AL, Lanford RE. **Cell surface assembly of lipoprotein(a) in primary cultures of baboon hepatocytes.** *J Biol Chem.* 1994;269(46):28716-23. — evidence that assembly takes place at the cell surface / extracellularly. The reason the model makes free apo(a) a separate state variable. [PMID 7961824](https://pubmed.ncbi.nlm.nih.gov/7961824/)
18. White AL, Hixson JE, Rainwater DL, Lanford RE. **Molecular basis for "null" lipoprotein(a) phenotypes and the influence of apolipoprotein(a) size on plasma lipoprotein(a) level in the baboon.** *J Biol Chem.* 1994;269(12):9060-6. — presecretory degradation of large isoforms. The direct basis of the `KDEGP` term in the model. [PMID 8132644](https://pubmed.ncbi.nlm.nih.gov/8132644/)
19. Brunner C, Lobentanz EM, Pethö-Schramm A, et al. **The number of identical kringle IV repeats in apolipoprotein(a) affects its processing and secretion by HepG2 cells.** *J Biol Chem.* 1996;271(50):32403-10. — copy number → ER residence time → secretion efficiency. The basis of the model's `SECEFF = KSZ^3/(KSZ^3+n^3)`. [PMID 8943305](https://pubmed.ncbi.nlm.nih.gov/8943305/)
20. Koschinsky ML, Côté GP, Gabel B, van der Hoek YY. **Identification of the cysteine residue in apolipoprotein(a) that mediates extracellular coupling with apolipoprotein B-100.** *J Biol Chem.* 1993;268(26):19819-25. — Cys4057. STEP 2 in the model. [PMID 8366118](https://pubmed.ncbi.nlm.nih.gov/8366118/)
21. Gabel BR, Koschinsky ML. **Sequences within apolipoprotein(a) kringle IV types 6-8 bind directly to low-density lipoprotein and mediate noncovalent association of apolipoprotein(a) with apolipoprotein B-100.** *Biochemistry.* 1998;37(21):7892-8. — STEP 1 (non-covalent docking). The muvalaplin target. [PMID 9601053](https://pubmed.ncbi.nlm.nih.gov/9601053/)
22. Becker L, McLeod RS, Marcovina SM, Yao Z, Koschinsky ML. **Identification of a critical lysine residue in apolipoprotein B-100 that mediates noncovalent interaction with apolipoprotein(a).** *J Biol Chem.* 2001;276(39):36155-62. [PMID 11469568](https://pubmed.ncbi.nlm.nih.gov/11469568/)
23. Youssef A, Clark JR, Koschinsky ML, Boffa MB. **Lipoprotein(a): expanding our knowledge of aortic valve narrowing.** *Trends Cardiovasc Med.* 2021;31(5):305-11. [PMID 32565142](https://pubmed.ncbi.nlm.nih.gov/32565142/)

## 4. Kinetics — production, not clearance, is what determines it

24. Rader DJ, Cain W, Zech LA, Usher D, Brewer HB Jr. **Variation in lipoprotein(a) concentrations among individuals with the same apolipoprotein(a) isoform is determined by the rate of lipoprotein(a) production.** *J Clin Invest.* 1993;91(2):443-7. — the structural premise of this whole model. [PMID 8432853](https://pubmed.ncbi.nlm.nih.gov/8432853/)
25. Rader DJ, Cain W, Ikewaki K, et al. **The inverse association of plasma lipoprotein(a) concentrations with apolipoprotein(a) isoform size is not due to differences in Lp(a) catabolism but to differences in production rate.** *J Clin Invest.* 1994;93(6):2758-63. — the FCR is independent of isoform size. [PMID 8201014](https://pubmed.ncbi.nlm.nih.gov/8201014/)
26. Frischmann ME, Ikewaki K, Trenkwalder E, et al. **In vivo stable-isotope kinetic study suggests intracellular assembly of lipoprotein(a).** *Atherosclerosis.* 2012;225(2):322-7. — contrary evidence on the site of assembly. The model adopts extracellular assembly but states this uncertainty explicitly. [PMID 23099120](https://pubmed.ncbi.nlm.nih.gov/23099120/)
27. Croyal M, Blanchard V, Ouguerram K, et al. **VLDL (very-low-density lipoprotein)-apo E (apolipoprotein E) may influence Lp(a) (lipoprotein [a]) synthesis or assembly.** *Arterioscler Thromb Vasc Biol.* 2020;40(3):819-29. [PMID 31941383](https://pubmed.ncbi.nlm.nih.gov/31941383/)
28. Chan DC, Watts GF, Coll B, Wasserman SM, Marcovina SM, Barrett PHR. **Lipoprotein(a) particle production as a determinant of plasma lipoprotein(a) concentration across varying apolipoprotein(a) isoform sizes and background cholesterol-lowering therapy.** *J Am Heart Assoc.* 2019;8(7):e011781. [PMID 30897997](https://pubmed.ncbi.nlm.nih.gov/30897997/)
29. Watts GF, Chan DC, Somaratne R, et al. **Controlled study of the effect of proprotein convertase subtilisin-kexin type 9 inhibition with evolocumab on lipoprotein(a) particle kinetics.** *Eur Heart J.* 2018;39(27):2577-85. — shows directly that the mechanism by which evolocumab lowers Lp(a) is increased catabolism. The basis of the `KLDLR_LPA` term in the model. [PMID 29566128](https://pubmed.ncbi.nlm.nih.gov/29566128/)
30. Reyes-Soffer G, Ginsberg HN, Ramakrishnan R. **The metabolism of lipoprotein(a): a story of dogma, extrapolation and controversy.** *Curr Opin Lipidol.* 2017;28(1):11-15. [PMID 27906712](https://pubmed.ncbi.nlm.nih.gov/27906712/)

## 5. Catabolic routes · receptors

31. Cain WJ, Millar JS, Himebauch AS, et al. **Lipoprotein(a) is cleared from the plasma primarily by the liver in a process mediated by apolipoprotein(a).** *J Lipid Res.* 2005;46(12):2681-91. [PMID 16150825](https://pubmed.ncbi.nlm.nih.gov/16150825/)
32. Romagnuolo R, Scipione CA, Boffa MB, Marcovina SM, Seidah NG, Koschinsky ML. **Lipoprotein(a) catabolism is regulated by proprotein convertase subtilisin/kexin type 9 through the low density lipoprotein receptor.** *J Biol Chem.* 2015;290(18):11649-62. — the existence of LDLR-mediated catabolism, and its limits. [PMID 25778403](https://pubmed.ncbi.nlm.nih.gov/25778403/)
33. Yang XP, Amar MJ, Vaisman B, et al. **Scavenger receptor-BI is a receptor for lipoprotein(a).** *J Lipid Res.* 2013;54(9):2450-7. [PMID 23812625](https://pubmed.ncbi.nlm.nih.gov/23812625/)
34. Sharma M, Redpath GM, Williams MJ, McCormick SP. **Recycling of apolipoprotein(a) after PlgRKT-mediated endocytosis of lipoprotein(a).** *Circ Res.* 2017;120(7):1091-1102. [PMID 28003219](https://pubmed.ncbi.nlm.nih.gov/28003219/)
35. Nielsen MB, Çolak Y, Benn M, Nordestgaard BG. **Plasma lipoprotein(a) and risk of aortic valve stenosis: the Copenhagen General Population Study.** *Eur Heart J.* 2019;40(38):3148-56. [PMID 31306481](https://pubmed.ncbi.nlm.nih.gov/31306481/)
36. Nioi P, Sigurdsson A, Thorleifsson G, et al. **Variant ASGR1 associated with a reduced risk of coronary artery disease.** *N Engl J Med.* 2016;374(22):2131-41. — loss of ASGR1 function lowers CAD risk. It is also the entry route for GalNAc drugs. [PMID 27192541](https://pubmed.ncbi.nlm.nih.gov/27192541/)
37. Kronenberg F, Utermann G, Dieplinger H. **Lipoprotein(a) in renal disease.** *Am J Kidney Dis.* 1996;27(1):1-25. — the renal catabolic route. The `RENF` term in the model. [PMID 8546112](https://pubmed.ncbi.nlm.nih.gov/8546112/)
38. Kostner KM, Maurer G, Huber K, et al. **Urinary excretion of apo(a) fragments. Role in apo(a) catabolism.** *Arterioscler Thromb Vasc Biol.* 1996;16(8):905-11. — the `FRAG` compartment in the model. [PMID 8696954](https://pubmed.ncbi.nlm.nih.gov/8696954/)

## 6. Measurement · isoform bias · units

39. Marcovina SM, Albers JJ, Gabel B, Koschinsky ML, Gaur VP. **Effect of the number of apolipoprotein(a) kringle 4 domains on immunochemical measurements of lipoprotein(a).** *Clin Chem.* 1995;41(2):246-55. — **the key evidence for the §5 (measurement) cluster of this model.** A polyclonal-antibody assay calibrated on a single isoform over-reports large isoforms and under-reports small ones. The model's `EPIT = (n+10)/(n_cal+10)`. [PMID 7533064](https://pubmed.ncbi.nlm.nih.gov/7533064/)
40. Marcovina SM, Albers JJ, Scanu AM, et al. **Use of a reference material proposed by the International Federation of Clinical Chemistry and Laboratory Medicine to evaluate analytical methods for the determination of plasma lipoprotein(a).** *Clin Chem.* 2000;46(12):1956-67. [PMID 11106328](https://pubmed.ncbi.nlm.nih.gov/11106328/)
41. Marcovina SM, Albers JJ. **Lipoprotein (a) measurements for clinical application.** *J Lipid Res.* 2016;57(4):526-37. — the basis of the recommendation to report in nmol/L. [PMID 26637279](https://pubmed.ncbi.nlm.nih.gov/26637279/)
42. Tsimikas S, Fazio S, Viney NJ, Xia S, Witztum JL, Marcovina SM. **Relationship of lipoprotein(a) molar concentration and mass according to lipoprotein(a) thresholds and apolipoprotein(a) isoform size.** *J Clin Lipidol.* 2018;12(5):1313-23. — the paper showing by measurement that no single conversion factor exists between mg/dL and nmol/L. Corresponds directly to RUN 3 of the model. [PMID 30037539](https://pubmed.ncbi.nlm.nih.gov/30037539/)
43. Kronenberg F, Mora S, Stroes ESG, et al. **Lipoprotein(a) in atherosclerotic cardiovascular disease and aortic stenosis: a European Atherosclerosis Society consensus statement.** *Eur Heart J.* 2022;43(39):3925-46. — the source of the ≥125 nmol/L / ≥50 mg/dL thresholds. [PMID 36036785](https://pubmed.ncbi.nlm.nih.gov/36036785/)
44. Koschinsky ML, Bajaj A, Boffa MB, et al. **A focused update to the 2019 NLA scientific statement on use of lipoprotein(a) in clinical practice.** *J Clin Lipidol.* 2024;18(3):e308-19. [PMID 38565461](https://pubmed.ncbi.nlm.nih.gov/38565461/)
45. Dahlén GH. **Incidence of Lp(a) lipoprotein among populations.** In: Scanu AM, ed. *Lipoprotein(a).* Academic Press; 1990:151-73. — the original source of the 0.30 coefficient used in the LDL-C correction.
46. Yeang C, Witztum JL, Tsimikas S. **'LDL-C' = LDL-C + Lp(a)-C: implications of achieved ultra-low LDL-C levels in the proprotein convertase subtilisin/kexin type 9 era of potent LDL-C lowering.** *Curr Opin Lipidol.* 2015;26(3):169-78. — the basis of the model's `LDLC_MEAS = LDLC_TRUE + LPAC`. [PMID 25943842](https://pubmed.ncbi.nlm.nih.gov/25943842/)
47. Yeang C, Witztum JL, Tsimikas S. **Novel method for quantification of lipoprotein(a)-cholesterol: implications for improving accuracy of LDL-C measurements.** *J Lipid Res.* 2021;62:100053. — evidence that the cholesterol fraction of Lp(a) mass may be ~17-25% rather than 30%. The uncertainty in the model's `FCHOL`. [PMID 33636162](https://pubmed.ncbi.nlm.nih.gov/33636162/)
48. Ruhaak LR, Cobbaert CM. **Quantifying apolipoprotein(a) in the era of proteoforms and precision medicine.** *Clin Chim Acta.* 2020;511:260-8. [PMID 33096035](https://pubmed.ncbi.nlm.nih.gov/33096035/)

## 7. Pathophysiology — the three effector arms

### 7a. Atherogenesis · matrix retention

49. Nielsen LB, Stender S, Kjeldsen K, Nordestgaard BG. **Specific accumulation of lipoprotein(a) in balloon-injured rabbit aorta in vivo.** *Circ Res.* 1996;78(4):615-26. — selective intimal accumulation of Lp(a). The `AVID` parameter in the model. [PMID 8635219](https://pubmed.ncbi.nlm.nih.gov/8635219/)
50. van der Hoek YY, Sangrar W, Côté GP, Kastelein JJ, Koschinsky ML. **Binding of recombinant apolipoprotein(a) to extracellular matrix proteins.** *Arterioscler Thromb.* 1994;14(11):1792-8. — fibronectin binding. [PMID 7947605](https://pubmed.ncbi.nlm.nih.gov/7947605/)
51. Boffa MB, Koschinsky ML. **Oxidized phospholipids as a unifying theory for lipoprotein(a) and cardiovascular disease.** *Nat Rev Cardiol.* 2019;16(5):305-18. [PMID 30675027](https://pubmed.ncbi.nlm.nih.gov/30675027/)
52. Björnson E, Adiels M, Taskinen MR, et al. **Lipoprotein(a) is markedly more atherogenic than LDL: an apolipoprotein B-based genetic analysis.** *J Am Coll Cardiol.* 2024;83(3):385-95. — the per-particle risk is far greater than that of LDL. The quantitative basis of the model's `AVID` · `WOX`. [PMID 38199713](https://pubmed.ncbi.nlm.nih.gov/38199713/)

### 7b. Oxidised phospholipids · inflammation

53. Bergmark C, Dewan A, Orsoni A, et al. **A novel function of lipoprotein(a) as a preferential carrier of oxidized phospholipids in human plasma.** *J Lipid Res.* 2008;49(10):2230-9. — the basis of ARM 2. [PMID 18594118](https://pubmed.ncbi.nlm.nih.gov/18594118/)
54. Tsimikas S, Brilakis ES, Miller ER, et al. **Oxidized phospholipids, Lp(a) lipoprotein, and coronary artery disease.** *N Engl J Med.* 2005;353(1):46-57. [PMID 16000355](https://pubmed.ncbi.nlm.nih.gov/16000355/)
55. Leibundgut G, Scipione C, Yin H, et al. **Determinants of binding of oxidized phospholipids on apolipoprotein(a) and lipoprotein(a).** *J Lipid Res.* 2013;54(10):2815-30. — OxPL is covalently bound to KIV-10. One site per particle. [PMID 23948545](https://pubmed.ncbi.nlm.nih.gov/23948545/)
56. van der Valk FM, Bekkering S, Kroon J, et al. **Oxidized phospholipids on lipoprotein(a) elicit arterial wall inflammation and an inflammatory monocyte response in humans.** *Circulation.* 2016;134(8):611-24. — trained immunity. Reversible with an apo(a) ASO. The `MONO` state variable in the model. [PMID 27496857](https://pubmed.ncbi.nlm.nih.gov/27496857/)
57. Müller N, Schulte DM, Türk K, et al. **IL-6 blockade by monoclonal antibodies inhibits apolipoprotein (a) expression and lipoprotein (a) synthesis in humans.** *J Lipid Res.* 2015;56(5):1034-42. — the IL-6 response element of the LPA promoter. Lp(a) falls under tocilizumab. The direct basis of the `FIL6` term in the model. [PMID 25713100](https://pubmed.ncbi.nlm.nih.gov/25713100/)
58. Wade DP, Clarke JG, Lindahl GE, et al. **5' control regions of the apolipoprotein(a) gene and members of the related plasminogen gene family.** *Proc Natl Acad Sci USA.* 1993;90(4):1369-73. — promoter structure. [PMID 8433995](https://pubmed.ncbi.nlm.nih.gov/8433995/)
59. Ridker PM, Devalaraja M, Baeres FMM, et al. **IL-6 inhibition with ziltivekimab in patients at high atherosclerotic risk (RESCUE): a double-blind, randomised, placebo-controlled, phase 2 trial.** *Lancet.* 2021;397(10289):2060-9. — hsCRP -92%, Lp(a) reduced. The data corresponding to scenario 15 of the model. [PMID 34015342](https://pubmed.ncbi.nlm.nih.gov/34015342/)
60. Ridker PM, Everett BM, Thuren T, et al. **Antiinflammatory therapy with canakinumab for atherosclerotic disease (CANTOS).** *N Engl J Med.* 2017;377(12):1119-31. — MACE reduced with no change in lipids. The reason the model's 'vulnerability' term has to exist independently. [PMID 28845751](https://pubmed.ncbi.nlm.nih.gov/28845751/)

### 7c. Antifibrinolysis — including the limits of the evidence

61. Miles LA, Fless GM, Levin EG, Scanu AM, Plow EF. **A potential basis for the thrombotic risks associated with lipoprotein(a).** *Nature.* 1989;339(6222):301-3. [PMID 2542796](https://pubmed.ncbi.nlm.nih.gov/2542796/)
62. Hancock MA, Boffa MB, Marcovina SM, Nesheim ME, Koschinsky ML. **Inhibition of plasminogen activation by lipoprotein(a): critical domains in apolipoprotein(a) and mechanism of inhibition on fibrin and degraded fibrin surfaces.** *J Biol Chem.* 2003;278(26):23260-9. — the strong lysine-binding site of KIV-10. One per particle. [PMID 12697750](https://pubmed.ncbi.nlm.nih.gov/12697750/)
63. Boffa MB. **Beyond fibrinolysis: the confounding role of Lp(a) in thrombosis.** *Atherosclerosis.* 2022;349:72-81. [PMID 35606082](https://pubmed.ncbi.nlm.nih.gov/35606082/)
64. Kamstrup PR, Tybjærg-Hansen A, Nordestgaard BG. **Genetic evidence that lipoprotein(a) associates with atherosclerotic stenosis rather than venous thrombosis.** *Arterioscler Thromb Vasc Biol.* 2012;32(7):1732-41. — **the reason the model sets the weight of ARM 3 (WARM3) low.** However strong the in vitro biology, there is no causal association with venous thrombosis. [PMID 22516069](https://pubmed.ncbi.nlm.nih.gov/22516069/)
65. Helgadottir A, Gretarsdottir S, Thorleifsson G, et al. **Apolipoprotein(a) genetic sequence variants associated with systemic atherosclerosis and coronary atherosclerotic burden but not with venous thromboembolism.** *J Am Coll Cardiol.* 2012;60(8):722-9. [PMID 22898070](https://pubmed.ncbi.nlm.nih.gov/22898070/)

### 7d. Calcific aortic valve disease

66. Thanassoulis G, Campbell CY, Owens DS, et al. **Genetic associations with valvular calcification and aortic stenosis.** *N Engl J Med.* 2013;368(6):503-12. — the causal association of LPA rs10455872 with aortic valve calcification. [PMID 23388002](https://pubmed.ncbi.nlm.nih.gov/23388002/)
67. Bouchareb R, Mahmut A, Nsaibia MJ, et al. **Autotaxin derived from lipoprotein(a) and valve interstitial cells promotes inflammation and mineralization of the aortic valve.** *Circulation.* 2015;132(8):677-90. — autotaxin → LysoPA → osteogenic transition. The skeleton of the valve module in the model. [PMID 26224810](https://pubmed.ncbi.nlm.nih.gov/26224810/)
68. Zheng KH, Tsimikas S, Pawade T, et al. **Lipoprotein(a) and oxidized phospholipids promote valve calcification in patients with aortic stenosis.** *J Am Coll Cardiol.* 2019;73(17):2150-62. — microcalcification progression confirmed by ¹⁸F-NaF PET. [PMID 31047003](https://pubmed.ncbi.nlm.nih.gov/31047003/)
69. Cowell SJ, Newby DE, Prescott RJ, et al. **A randomized trial of intensive lipid-lowering therapy in calcific aortic stenosis (SALTIRE).** *N Engl J Med.* 2005;352(23):2389-97. — **a negative result.** The phenomenon the model's self-perpetuation term (KSELFR) reproduces. [PMID 15944423](https://pubmed.ncbi.nlm.nih.gov/15944423/)
70. Rossebø AB, Pedersen TR, Boman K, et al. **Intensive lipid lowering with simvastatin and ezetimibe in aortic stenosis (SEAS).** *N Engl J Med.* 2008;359(13):1343-56. — a negative result. [PMID 18765433](https://pubmed.ncbi.nlm.nih.gov/18765433/)
71. Chan KL, Teo K, Dumesnil JG, Ni A, Tam J; ASTRONOMER Investigators. **Effect of lipid lowering with rosuvastatin on progression of aortic stenosis (ASTRONOMER).** *Circulation.* 2010;121(2):306-14. — a negative result. [PMID 20048204](https://pubmed.ncbi.nlm.nih.gov/20048204/)

## 8. The Lp(a) effect of established lipid drugs

72. Tsimikas S, Gordts PLSM, Nora C, Yeang C, Witztum JL. **Statin therapy increases lipoprotein(a) levels.** *Eur Heart J.* 2020;41(24):2275-84. — meta-analysis. What scenario 4 of the model has to reproduce. [PMID 31111151](https://pubmed.ncbi.nlm.nih.gov/31111151/)
73. de Boer LM, Oorthuys AOJ, Wiegman A, et al. **Statin therapy and lipoprotein(a) levels: a systematic review and meta-analysis.** *Eur J Prev Cardiol.* 2022;29(5):779-92. [PMID 34849691](https://pubmed.ncbi.nlm.nih.gov/34849691/)
74. O'Donoghue ML, Fazio S, Giugliano RP, et al. **Lipoprotein(a), PCSK9 inhibition, and cardiovascular risk: insights from the FOURIER trial.** *Circulation.* 2019;139(12):1483-92. [PMID 30586750](https://pubmed.ncbi.nlm.nih.gov/30586750/)
75. Bittner VA, Szarek M, Aylward PE, et al. **Effect of alirocumab on lipoprotein(a) and cardiovascular risk after acute coronary syndrome (ODYSSEY OUTCOMES).** *J Am Coll Cardiol.* 2020;75(2):133-44. — the fall in Lp(a) contributes to the MACE reduction independently of LDL-C. [PMID 31948641](https://pubmed.ncbi.nlm.nih.gov/31948641/)
76. Ray KK, Wright RS, Kallend D, et al. **Two phase 3 trials of inclisiran in patients with elevated LDL cholesterol (ORION-10/-11).** *N Engl J Med.* 2020;382(16):1507-19. — Lp(a) about -19 ~ -26%. [PMID 32187462](https://pubmed.ncbi.nlm.nih.gov/32187462/)
77. Awad K, Mikhailidis DP, Katsiki N, Muntner P, Banach M; Lipid and Blood Pressure Meta-Analysis Collaboration Group. **Effect of ezetimibe monotherapy on plasma lipoprotein(a) concentrations in patients with primary hypercholesterolemia: a systematic review and meta-analysis.** *Drugs.* 2018;78(4):453-62. — the small (about -7%) Lp(a) reduction with ezetimibe. **A prediction the model reproduced without post hoc adjustment** (see §13). [PMID 29396832](https://pubmed.ncbi.nlm.nih.gov/29396832/)
78. Boden WE, Probstfield JL, Anderson T, et al; AIM-HIGH Investigators. **Niacin in patients with low HDL cholesterol levels receiving intensive statin therapy.** *N Engl J Med.* 2011;365(24):2255-67. — a negative result. [PMID 22085343](https://pubmed.ncbi.nlm.nih.gov/22085343/)
79. HPS2-THRIVE Collaborative Group. **Effects of extended-release niacin with laropiprant in high-risk patients.** *N Engl J Med.* 2014;371(3):203-12. — a negative result. The model explains it not as 'insufficient potency' but as 'insufficient absolute reduction'. [PMID 25014686](https://pubmed.ncbi.nlm.nih.gov/25014686/)
80. Nicholls SJ, Ditmarsch M, Kastelein JJ, et al. **Lipid lowering effects of the CETP inhibitor obicetrapib in combination with high-intensity statins (ROSE2).** *Nat Med.* 2022;28(8):1672-8. [PMID 35953719](https://pubmed.ncbi.nlm.nih.gov/35953719/)
81. Nicholls SJ, Nelson AJ, Ditmarsch M, et al. **Obicetrapib on top of maximally tolerated lipid-modifying therapies in participants with or at high risk for ASCVD (BROADWAY).** *Nat Med.* 2025;31(2):500-8. [PMID 39653774](https://pubmed.ncbi.nlm.nih.gov/39653774/)
82. Roeseler E, Julius U, Heigl F, et al; Pro(a)LiFe-Study Group. **Lipoprotein apheresis for lipoprotein(a)-associated cardiovascular disease: prospective 5 years of follow-up and apo(a) characterization.** *Arterioscler Thromb Vasc Biol.* 2016;36(9):2019-27. [PMID 27417585](https://pubmed.ncbi.nlm.nih.gov/27417585/)
83. Waldmann E, Parhofer KG. **Lipoprotein apheresis to treat elevated lipoprotein(a).** *J Lipid Res.* 2016;57(10):1751-7. — acute removal of 60-70%, interval mean 30-35%. Scenario 16 of the model. [PMID 26658193](https://pubmed.ncbi.nlm.nih.gov/26658193/)

## 9. RNA-directed therapeutics

84. Tsimikas S, Viney NJ, Hughes SG, et al. **Antisense therapy targeting apolipoprotein(a): a randomised, double-blind, placebo-controlled phase 1 study.** *Lancet.* 2015;386(10002):1472-83. [PMID 26210642](https://pubmed.ncbi.nlm.nih.gov/26210642/)
85. Viney NJ, van Capelleveen JC, Geary RS, et al. **Antisense oligonucleotides targeting apolipoprotein(a) in people with raised lipoprotein(a): two randomised, double-blind, placebo-controlled, dose-ranging trials.** *Lancet.* 2016;388(10057):2239-53. [PMID 27665230](https://pubmed.ncbi.nlm.nih.gov/27665230/)
86. Tsimikas S, Karwatowska-Prokopczuk E, Gouni-Berthold I, et al. **Lipoprotein(a) reduction in persons with cardiovascular disease (pelacarsen phase 2b).** *N Engl J Med.* 2020;382(3):244-55. — 80 mg once monthly → about -80%. The target of scenario 8 of the model. [PMID 31893580](https://pubmed.ncbi.nlm.nih.gov/31893580/)
87. O'Donoghue ML, Rosenson RS, Gencer B, et al; OCEAN(a)-DOSE Trial Investigators. **Small interfering RNA to reduce lipoprotein(a) in cardiovascular disease.** *N Engl J Med.* 2022;387(20):1855-64. — olpasiran 75 mg q12w → about -95 ~ -101% placebo-adjusted. [PMID 36342163](https://pubmed.ncbi.nlm.nih.gov/36342163/)
88. Nissen SE, Wolski K, Watts GF, et al. **Single ascending and multiple-dose trial of zerlasiran, a short interfering RNA targeting lipoprotein(a).** *JAMA.* 2024;331(17):1472-81. [PMID 38583084](https://pubmed.ncbi.nlm.nih.gov/38583084/)
89. Nissen SE, Linnebjerg H, Shen X, et al. **Lepodisiran, an extended-duration short interfering RNA targeting lipoprotein(a): a randomized dose-ascending clinical trial.** *JAMA.* 2023;330(21):2075-83. [PMID 37952254](https://pubmed.ncbi.nlm.nih.gov/37952254/)
90. Nissen SE, Wang Q, Nicholls SJ, et al; ALPACA Investigators. **Lepodisiran in adults with elevated lipoprotein(a).** *N Engl J Med.* 2025;392(19):1892-1902. — a single 608 mg dose gives >90% reduction for about 48 weeks. [PMID 40162914](https://pubmed.ncbi.nlm.nih.gov/40162914/)
91. Nissen SE, Linnebjerg H, Shen X, et al. **Muvalaplin, an oral small molecule inhibitor of lipoprotein(a) formation: a randomized clinical trial (phase 1).** *JAMA.* 2023;330(11):1042-53. — assembly inhibition as a new target. A rise in free apo(a) was observed. [PMID 37638695](https://pubmed.ncbi.nlm.nih.gov/37638695/)
92. Nicholls SJ, Ni W, Rhodes GM, et al. **Oral muvalaplin for lowering of lipoprotein(a): a randomized clinical trial (KRAKEN).** *JAMA.* 2025;333(3):222-31. — **-85.8% by the intact-Lp(a) assay, a smaller reduction by the traditional apo(a) assay.** What §C of the model (the muvalaplin assay discrepancy) addresses. [PMID 39556753](https://pubmed.ncbi.nlm.nih.gov/39556753/)
93. Crooke ST, Baker BF, Crooke RM, Liang XH. **Antisense technology: an overview and prospectus.** *Nat Rev Drug Discov.* 2021;20(6):427-53. — GalNAc-ASGR1 hepatocyte targeting. The PK structure of the model. [PMID 33762737](https://pubmed.ncbi.nlm.nih.gov/33762737/)
94. Springer AD, Dowdy SF. **GalNAc-siRNA conjugates: leading the way for delivery of RNAi therapeutics.** *Nucleic Acid Ther.* 2018;28(3):109-18. — the long half-life of RISC loading. The basis of the `SIR_RISC` compartment in the model. [PMID 29792572](https://pubmed.ncbi.nlm.nih.gov/29792572/)

## 10. Epidemiology · risk quantification

95. Erqou S, Kaptoge S, Perry PL, et al; Emerging Risk Factors Collaboration. **Lipoprotein(a) concentration and the risk of coronary heart disease, stroke, and nonvascular mortality.** *JAMA.* 2009;302(4):412-23. — HR 1.13 per 1 SD. [PMID 19622820](https://pubmed.ncbi.nlm.nih.gov/19622820/)
96. Burgess S, Ference BA, Staley JR, et al. **Association of LPA variants with risk of coronary disease and the implications for lipoprotein(a)-lowering therapies: a Mendelian randomization analysis.** *JAMA Cardiol.* 2018;3(7):619-27. — **a lifelong reduction of 101.5 mg/dL ≈ a 38.67 mg/dL reduction in LDL-C.** The source of the risk-translation coefficient in the model. [PMID 29926099](https://pubmed.ncbi.nlm.nih.gov/29926099/)
97. Lamina C, Kronenberg F; Lp(a)-GWAS-Consortium. **Estimation of the required lipoprotein(a)-lowering therapeutic effect size for reduction in coronary heart disease outcomes: a Mendelian randomization analysis.** *JAMA Cardiol.* 2019;4(6):575-9. — the absolute reduction required. The argument behind RUN 4 of the model. [PMID 31017644](https://pubmed.ncbi.nlm.nih.gov/31017644/)
98. Ference BA, Ginsberg HN, Graham I, et al. **Low-density lipoproteins cause atherosclerotic cardiovascular disease. 1. Evidence from genetic, epidemiologic, and clinical studies.** *Eur Heart J.* 2017;38(32):2459-72. — the difference in effect size between lifelong exposure and a short-term trial (about threefold). The basis of the two-time-constant structure of the model. [PMID 28444290](https://pubmed.ncbi.nlm.nih.gov/28444290/)
99. Cholesterol Treatment Trialists' (CTT) Collaboration. **Efficacy and safety of more intensive lowering of LDL cholesterol: a meta-analysis of data from 170,000 participants in 26 randomised trials.** *Lancet.* 2010;376(9753):1670-81. — the per-unit effect in a five-year trial. [PMID 21067804](https://pubmed.ncbi.nlm.nih.gov/21067804/)
100. Nordestgaard BG, Chapman MJ, Ray K, et al. **Lipoprotein(a) as a cardiovascular risk factor: current status.** *Eur Heart J.* 2010;31(23):2844-53. [PMID 20965889](https://pubmed.ncbi.nlm.nih.gov/20965889/)
101. Willeit P, Ridker PM, Nestel PJ, et al. **Baseline and on-statin treatment lipoprotein(a) levels for prediction of cardiovascular events: individual patient-data meta-analysis of statin outcome trials.** *Lancet.* 2018;392(10155):1311-20. — residual risk. [PMID 30293769](https://pubmed.ncbi.nlm.nih.gov/30293769/)
102. Berman AN, Biery DW, Besser SA, et al. **Lipoprotein(a) and major adverse cardiovascular events in patients with or without baseline atherosclerotic cardiovascular disease.** *J Am Coll Cardiol.* 2024;83(9):873-86. [PMID 38418000](https://pubmed.ncbi.nlm.nih.gov/38418000/)

## 11. Secondary causes · modifiers

103. Wanner C, Rader D, Bartens W, et al. **Elevated plasma lipoprotein(a) in patients with the nephrotic syndrome.** *Ann Intern Med.* 1993;119(4):263-9. — the model's `FNEPH`. [PMID 8328734](https://pubmed.ncbi.nlm.nih.gov/8328734/)
104. de Bruin TW, van Barlingen H, van Linde-Sibenius Trip M, van Vuurst de Vries AR, Akveld MJ, Erkelens DW. **Lipoprotein(a) and apolipoprotein B plasma concentrations in hypothyroid, euthyroid, and hyperthyroid subjects.** *J Clin Endocrinol Metab.* 1993;76(1):121-6. — the model's `FTHY`. [PMID 8421075](https://pubmed.ncbi.nlm.nih.gov/8421075/)
105. Anagnostis P, Galanis P, Chatzistergiou V, et al. **The effect of hormone replacement therapy and tibolone on lipoprotein(a) concentrations in postmenopausal women: a systematic review and meta-analysis.** *Maturitas.* 2017;99:27-36. — the model's `FSEX`. [PMID 28364865](https://pubmed.ncbi.nlm.nih.gov/28364865/)
106. Enkhmaa B, Anuurad E, Berglund L. **Lipoprotein (a): impact by ethnicity and environmental and medical conditions.** *J Lipid Res.* 2016;57(7):1111-25. [PMID 26637279](https://pubmed.ncbi.nlm.nih.gov/26637279/)
107. Burgess S, Davey Smith G. **Mendelian randomization implicates adiposity-related traits in the pathogenesis of lipoprotein(a) elevation.** *Eur J Prev Cardiol.* 2018;25(6):617-9. [PMID 29517303](https://pubmed.ncbi.nlm.nih.gov/29517303/)
108. Norata GD, Ballantyne CM, Catapano AL. **New therapeutic principles in dyslipidaemia: focus on LDL and Lp(a) lowering drugs.** *Eur Heart J.* 2013;34(24):1783-9. [PMID 23509225](https://pubmed.ncbi.nlm.nih.gov/23509225/)

## 12. Guidelines · ongoing outcome trials

109. Mach F, Baigent C, Catapano AL, et al. **2019 ESC/EAS Guidelines for the management of dyslipidaemias.** *Eur Heart J.* 2020;41(1):111-88. — the recommendation to measure Lp(a) once in a lifetime. [PMID 31504418](https://pubmed.ncbi.nlm.nih.gov/31504418/)
110. Reyes-Soffer G, Ginsberg HN, Berglund L, et al; American Heart Association. **Lipoprotein(a): a genetically determined, causal, and prevalent risk factor for atherosclerotic cardiovascular disease — an AHA scientific statement.** *Arterioscler Thromb Vasc Biol.* 2022;42(1):e48-60. [PMID 34647487](https://pubmed.ncbi.nlm.nih.gov/34647487/)
111. Tsimikas S, Moriarty PM, Stroes ES. **Emerging RNA therapeutics to lower blood levels of Lp(a): JACC Focus Seminar 2/4.** *J Am Coll Cardiol.* 2021;77(12):1576-89. [PMID 33766265](https://pubmed.ncbi.nlm.nih.gov/33766265/)
112. Tsimikas S, Karwatowska-Prokopczuk E, Yeang C, et al. **Rationale and design of the Lp(a)HORIZON trial: assessing the effect of pelacarsen on major cardiovascular events in patients with CVD and elevated Lp(a).** *Am Heart J.* 2025;282:9-18. — enrolment criterion Lp(a) ≥70 mg/dL. Corresponds directly to the design logic of RUN 4 of the model. [PMID 39674305](https://pubmed.ncbi.nlm.nih.gov/39674305/)
113. Nissen SE, Wolski K, Balog C, et al. **OCEAN(a)-Outcomes: rationale and design of a trial of olpasiran in patients with elevated Lp(a) and established ASCVD.** *Am Heart J.* 2024;278:1-10. — enrolment criterion ≥200 nmol/L. [PMID 39179148](https://pubmed.ncbi.nlm.nih.gov/39179148/)

---

## 13. Equation ↔ source map

| Model element | Equation / parameter | Source |
|---|---|---|
| Production determines the variation | `KCAT0` fixed · `KTL` back-calculated | 24, 25, 8 |
| Isoform size → secretion efficiency | `SECEFF = KSZ³/(KSZ³+n³)` | 18, 19, 7 |
| Two-step assembly · free apo(a) | `APOA_FR`, `ASSEM` | 17, 20, 21, 22 |
| Assembly saturates in the LDL substrate | `LDL_P/(KMLDL+LDL_P)` | 77 (ezetimibe barely changes Lp(a)) |
| LDLR-dependent catabolism is the minor route | `KLDLR_LPA/KCAT0 = 0.24` | 29, 32, 74, 75 |
| The two opposing terms of a statin | `ESRE` vs `ESTA` | 72, 73 |
| IL-6 response element of the LPA promoter | `FIL6` | 57, 58, 59 |
| One OxPL site per particle | `OXPL` | 53, 54, 55 |
| Trained monocytes | `MONO` | 56 |
| Matrix-binding affinity 2.5-fold | `AVID` | 49, 50, 52 |
| One lysine-binding site per particle | `PLGOCC` | 62 |
| Low weight on the antifibrinolysis arm | `WARM3 = 0.25` | 64, 65 |
| Autotaxin → LysoPA → osteogenesis | `ATXV`→`LYSOPA`→`VIC_OST` | 67, 68 |
| Self-perpetuating valve calcification | `KSELFR` | 69, 70, 71 |
| Antibody bias of the mass assay | `EPIT = (n+10)/(n_cal+10)` | 39, 40, 42 |
| LDL-C contamination · the Dahlén correction | `LDLC_MEAS`, `LDLC_CORR` | 45, 46, 47 |
| Lifelong exposure versus a five-year trial | `BSLOW` (slow) + `BFAST` (fast) | 96, 97, 98, 99 |
| GalNAc-ASO / siRNA PK | `PEL_LIV`, `SIR_RISC` | 93, 94, 86, 87 |
| Assembly inhibitor | `FMUV` | 21, 91, 92 |

### Observations the model reproduced without post hoc adjustment (predictions, not fits)

| Observation | Literature value | Model value |
|---|---|---|
| Lp(a) effect of ezetimibe | about −7% (77) | **−6.6%** |
| PCSK9 inhibitor Lp(a) | −25 ~ −30% (74, 75) | **−30.1%** |
| Anti-IL-6 (non-inflammatory patients) | −16 ~ −25% (59) | **−26.6%** |
| Anti-IL-6 (rheumatoid arthritis) | −37% (57) | **−33.5%** |
| hsCRP (RESCUE) | −92% (59) | **−92.5%** |
| Pelacarsen apoB | about −13% (86) | **−11.3%** |

### What this model cannot answer — the open questions it exposes

1. **The size of the assay discrepancy in KRAKEN.** The model reproduces the direction (the intact assay reports the
   larger reduction) but not the size (−85.8% versus −70%). To close it, the traditional apo(a)
   assay would have to have a **molar response to free apo(a) about 3.8 times that to apo(a) inside a particle**
   (parameter `RFREE`). This is a measurable quantity and, as far as is known, has not yet been
   reported. (91, 92)
2. **The site of assembly.** The extracellular (17) versus intracellular (26) dispute is unresolved, and the model has
   adopted the extracellular. If intracellular assembly is right, the predicted rise in free apo(a) with muvalaplin changes.
3. **The cholesterol fraction of Lp(a) mass.** The classical 0.30 (45) versus the recent 0.17-0.25 (47).
   The model's `FCHOL` exposes this uncertainty as it stands, and whether the Dahlén correction over-corrects
   depends on it.
4. **Whether blocked assembly returns the spared apoB to plasma LDL.** The model's `FRECY`.
   The size of the apoB reduction in the pelacarsen trial (86) suggests `FRECY` is close to 0.
