# Acne vulgaris — QSP model references
### Reference list for `acn_qsp_model.dot` / `acn_mrgsolve_model.R` / `acn_shiny_app.R`

Every PMID in this list was looked up through the NCBI E-utilities and its title, journal and year confirmed.
The grey note under each entry records **which part of the model that paper supports**.
(The model parameters are approximations informed by the values reported in these papers; they have not been through a formal fitting.)


## A. Overview · epidemiology · disease burden

1. Williams HC, Dellavalle RP, Garner S. Acne vulgaris. *Lancet*. 2012;379:361-72. [PMID 21880356](https://pubmed.ncbi.nlm.nih.gov/21880356/)  
   <sub>A Lancet seminar. The standard summary of prevalence, natural history and treatment stratification. Used to define the model's reference 'moderate' phenotype.</sub>
2. Moradi Tuchayi S, Makrantonaki E, Ganceviciene R, et al. Acne vulgaris. *Nat Rev Dis Primers*. 2015;1:15033. [PMID 27227877](https://pubmed.ncbi.nlm.nih.gov/27227877/)  
   <sub>Nature Reviews Disease Primers. The structure of the four pathogenic pillars (excess sebum · hyperkeratinisation · C. acnes · inflammation) is itself the skeleton of this model.</sub>
3. Yuan R, Long H, Li T, et al. Comparison of the Burden of Acne Vulgaris in China vs Globally: Insights from the Global Burden of Disease Study (2021) and Projections to 2050. *Acta Derm Venereol*. 2026;106. [PMID 42051025](https://pubmed.ncbi.nlm.nih.gov/42051025/)  
   <sub>A GBD-based comparison of disease burden. The demographic background to the adolescent peak and the persistent adult form.</sub>
4. Samuels DV, Rosenthal R, Lin R, et al. Acne vulgaris and risk of depression and anxiety: A meta-analytic review. *J Am Acad Dermatol*. 2020;83:532-541. [PMID 32088269](https://pubmed.ncbi.nlm.nih.gov/32088269/)  
   <sub>Meta-analysis of the risk of depression and anxiety. The QoL context of the Shiny app and the argument for early treatment.</sub>
5. Veldi VDK, Metta AK, Metta S, et al. Living With Acne Vulgaris in Young Adults: A Holistic Examination of Its Impact on Quality of Life Using the Dermatology Life Quality Index (DLQI). *Cureus*. 2025;17:e77167. [PMID 39925569](https://pubmed.ncbi.nlm.nih.gov/39925569/)  
   <sub>A qualitative and quantitative review of the quality-of-life impact of acne in young people.</sub>
6. Salari Y, Latt M, Lau E, et al. The Prevalence and Burden of Truncal Acne. *Dermatol Ther (Heidelb)*. 2026;16:3245-3251. [PMID 41963698](https://pubmed.ncbi.nlm.nih.gov/41963698/)  
   <sub>Prevalence and burden of truncal acne. Background to the trifarotene truncal-indication scenario.</sub>
7. Friedman GD. Twin studies of disease heritability based on medical records: application to acne vulgaris. *Acta Genet Med Gemellol (Roma)*. 1984;33:487-95. [PMID 6241419](https://pubmed.ncbi.nlm.nih.gov/6241419/)  
   <sub>A medical-record-based twin study — the heritability of acne. Grounds for introducing the model's constitutional parameters SEVX/LESX.</sub>
8. Maxwell J, Mitchell BL, Du-Harpur X, et al. Genome-wide association meta-regression identifies stem cell lineage orchestration as a key driver of acne risk. *medRxiv*. 2025. [PMID 40666320](https://pubmed.ncbi.nlm.nih.gov/40666320/)  
   <sub>An acne GWAS meta-regression — hair follicle stem cell lineage genes. Cluster 1 of the mechanistic map.</sub>

## B. Pathogenesis in general · the microcomedone · keratinisation

9. Oulès B, Saurat JH. Strategic Targets in Acne, Update 2025: The Microcomedone Is Not Just a Plug, It Is an Egg. *Dermatology*. 2026;242:8-14. [PMID 40924653](https://pubmed.ncbi.nlm.nih.gov/40924653/)  
   <sub>A 2025 update. 'The microcomedone is not a simple plug' — the direct grounds for holding MC as a separate state variable in the model.</sub>
10. Jeremy AH, Holland DB, Roberts SG, et al. Inflammatory events are involved in acne lesion initiation. *J Invest Dermatol*. 2003;121:20-7. [PMID 12839559](https://pubmed.ncbi.nlm.nih.gov/12839559/)  
   <sub>Jeremy et al. The classic demonstration that inflammatory events already precede in a clinically normal follicle. Grounds for the de novo inflammatory papule term (KDENOVO).</sub>
11. Persson G, Johansson-Jänkänpää E, Ganceviciene R, et al. No evidence for follicular keratinocyte hyperproliferation in acne lesions as compared to autologous healthy hair follicles. *Exp Dermatol*. 2018;27:668-671. [PMID 29582469](https://pubmed.ncbi.nlm.nih.gov/29582469/)  
   <sub>The contrary study finding no evidence of follicular keratinocyte hyperproliferation in acne lesions. Why the model treats KER as an index of 'cohesion and differentiation abnormality' rather than of 'proliferation'.</sub>
12. Perisho K, Wertz PW, Madison KC, et al. Fatty acids of acylceramides from comedones and from the skin surface of acne patients and control subjects. *J Invest Dermatol*. 1988;90:350-3. [PMID 2964492](https://pubmed.ncbi.nlm.nih.gov/2964492/)  
   <sub>The fatty acid composition of acylceramides in comedones and on the acne skin surface — the measured grounds for the linoleic acid dilution hypothesis (state variable LA).</sub>
13. Das S, Reynolds RV. Recent advances in acne pathogenesis: implications for therapy. *Am J Clin Dermatol*. 2014;15:479-88. [PMID 25388823](https://pubmed.ncbi.nlm.nih.gov/25388823/)  
   <sub>Recent advances in the pathogenesis of acne and their therapeutic implications — an account of the interactions between the pillars.</sub>

## C. The androgen axis and sebaceous gland biology

14. Chen WC, Zouboulis CC. Hormones and the pilosebaceous unit. *Dermatoendocrinol*. 2009;1:81-6. [PMID 20224689](https://pubmed.ncbi.nlm.nih.gov/20224689/)  
   <sub>Hormones and the pilosebaceous unit. The general account of the ARS → SGM/LIP link.</sub>
15. Picardo M, Ottaviani M, Camera E, et al. Sebaceous gland lipids. *Dermatoendocrinol*. 2009;1:68-71. [PMID 20224686](https://pubmed.ncbi.nlm.nih.gov/20224686/)  
   <sub>A review of sebaceous gland lipids — the compositional proportions of squalene, wax esters and triglycerides (cluster 4 of the map).</sub>
16. Khondker L, Khan SI. Acne vulgaris related to androgens - a review. *Mymensingh Med J*. 2014;23:181-5. [PMID 24584396](https://pubmed.ncbi.nlm.nih.gov/24584396/)  
   <sub>A review of androgens and acne. Grounds for using the free androgen index (FAI).</sub>
17. Samson M, Labrie F, Zouboulis CC, et al. Biosynthesis of dihydrotestosterone by a pathway that does not require testosterone as an intermediate in the SZ95 sebaceous gland cell line. *J Invest Dermatol*. 2010;130:602-4. [PMID 19812596](https://pubmed.ncbi.nlm.nih.gov/19812596/)  
   <sub>The DHT biosynthetic route in sebocytes that does not pass through testosterone — the local steroidogenesis node.</sub>
18. Leyden J, Bergfeld W, Drake L, et al. A systemic type I 5 alpha-reductase inhibitor is ineffective in the treatment of acne vulgaris. *J Am Acad Dermatol*. 2004;50:443-7. [PMID 14988688](https://pubmed.ncbi.nlm.nih.gov/14988688/)  
   <sub>The trial in which a systemic type 1 5α-reductase inhibitor was ineffective in acne. Why the model leaves the E5ARI default at 0 and treats blockade of 5AR alone as a weak lever.</sub>
19. Cotterill JA, Cunliffe WJ, Williamson B. Severity of acne and sebum excretion rate. *Br J Dermatol*. 1971;85:93-4. [PMID 4254153](https://pubmed.ncbi.nlm.nih.gov/4254153/)  
   <sub>The relation between acne severity and the sebum excretion rate — calibration of the SER reference value (normal ~0.9 vs acne ~2 µg/cm²/min).</sub>
20. Pan J, Wang Q, Tu P. A Topical Medication of All-Trans Retinoic Acid Reduces Sebum Excretion Rate in Patients With Forehead Acne. *Am J Ther*. 2017;24:e207-e212. [PMID 26872139](https://pubmed.ncbi.nlm.nih.gov/26872139/)  
   <sub>The observed reduction in sebum excretion rate with topical all-trans retinoic acid — sets the upper bound on the size of the retinoid sebum effect.</sub>
21. Chiba K, Yoshizawa K, Makino I, et al. Comedogenicity of squalene monohydroperoxide in the skin after topical application. *J Toxicol Sci*. 2000;25:77-83. [PMID 10845185](https://pubmed.ncbi.nlm.nih.gov/10845185/)  
   <sub>The comedogenicity of squalene monohydroperoxide — the SQOX → comedone formation / oxidative colouring of the open comedone route.</sub>
22. Zouboulis CC, Oeff MK, Hiroi N, et al. Involvement of Pattern Recognition Receptors in the Direct Influence of Bacterial Components and Standard Antiacne Compounds on Human Sebaceous Gland Cells. *Skin Pharmacol Physiol*. 2021;34:19-29. [PMID 33601383](https://pubmed.ncbi.nlm.nih.gov/33601383/)  
   <sub>The involvement of pattern recognition receptors in sebocytes — the sebaceous gland is an immunologically active tissue, not merely a secretory organ.</sub>

## D. Nutrient-metabolic signalling (Insulin · IGF-1 · mTORC1 · FoxO1)

23. Cappel M, Mauger D, Thiboutot D. Correlation between serum levels of insulin-like growth factor 1, dehydroepiandrosterone sulfate, and dihydrotestosterone and acne lesion counts in adult women. *Arch Dermatol*. 2005;141:333-8. [PMID 15781674](https://pubmed.ncbi.nlm.nih.gov/15781674/)  
   <sub>The correlation of IGF-1, DHEAS and DHT with lesion count in adult women — sets the IGF1 → LIP weight (WIGF).</sub>
24. Smith R, Mann N, Mäkeläinen H, et al. A pilot study to determine the short-term effects of a low glycemic load diet on hormonal markers of acne: a nonrandomized, parallel, controlled feeding trial. *Mol Nutr Food Res*. 2008;52:718-26. [PMID 18496812](https://pubmed.ncbi.nlm.nih.gov/18496812/)  
   <sub>A randomised pilot study of a low glycaemic load diet — the upper bound on the effect size of the GLYLOAD parameter.</sub>
25. Adebamowo CA, Spiegelman D, Berkey CS, et al. Milk consumption and acne in teenaged boys. *J Am Acad Dermatol*. 2008;58:787-93. [PMID 18194824](https://pubmed.ncbi.nlm.nih.gov/18194824/)  
   <sub>Milk intake and acne in adolescent boys — the DAIRY parameter.</sub>
26. Juhl CR, Bergholdt HKM, Miller IM, et al. Dairy Intake and Acne Vulgaris: A Systematic Review and Meta-Analysis of 78,529 Children, Adolescents, and Young Adults. *Nutrients*. 2018;10. [PMID 30096883](https://pubmed.ncbi.nlm.nih.gov/30096883/)  
   <sub>Dairy intake and acne: a systematic review and meta-analysis in 78,000 subjects.</sub>
27. Melnik BC. Lifetime Impact of Cow's Milk on Overactivation of mTORC1: From Fetal to Childhood Overgrowth, Acne, Diabetes, Cancers, and Neurodegeneration. *Biomolecules*. 2021;11. [PMID 33803410](https://pubmed.ncbi.nlm.nih.gov/33803410/)  
   <sub>Overactivation of mTORC1 by milk — the leucine/BCAA node and the NUTR term.</sub>
28. Agamia NF, Hussein OM, Abdelmaksoud RE, et al. Effect of oral isotretinoin on the nucleo-cytoplasmic distribution of FoxO1 and FoxO3 proteins in sebaceous glands of patients with acne vulgaris. *Exp Dermatol*. 2018;27:1344-1351. [PMID 30240097](https://pubmed.ncbi.nlm.nih.gov/30240097/)  
   <sub>Human data showing that oral isotretinoin changes the nuclear-cytoplasmic distribution of FoxO1/FoxO3 — the grounds for the term by which isotretinoin directly inhibits LIP (mTORC1/SREBP) in the model (EISOLIP).</sub>

## E. Cutibacterium acnes — phylotypes · biofilm · metabolites

29. Fitz-Gibbon S, Tomida S, Chiu BH, et al. Propionibacterium acnes strain populations in the human skin microbiome associated with acne. *J Invest Dermatol*. 2013;133:2152-60. [PMID 23337890](https://pubmed.ncbi.nlm.nih.gov/23337890/)  
   <sub>The landmark study of Fitz-Gibbon et al. It is the 'strain composition', not the 'quantity' of the organism, that is associated with acne — why the model holds PHYLIA (the dominance of IA1) as a separate parameter.</sub>
30. McDowell A, Perry AL, Lambert PA, et al. A new phylogenetic group of Propionibacterium acnes. *J Med Microbiol*. 2008;57:218-224. [PMID 18201989](https://pubmed.ncbi.nlm.nih.gov/18201989/)  
   <sub>Description of a new phylogenetic group of C. acnes.</sub>
31. Dagnelie MA, Corvec S, Saint-Jean M, et al. Cutibacterium acnes phylotypes diversity loss: a trigger for skin inflammatory process. *J Eur Acad Dermatol Venereol*. 2019;33:2340-2348. [PMID 31299116](https://pubmed.ncbi.nlm.nih.gov/31299116/)  
   <sub>Direct evidence that the loss of phylotype diversity is the trigger of skin inflammation.</sub>
32. Dreno B, Dekio I, Baldwin H, et al. Acne microbiome: From phyla to phylotypes. *J Eur Acad Dermatol Venereol*. 2024;38:657-664. [PMID 37777343](https://pubmed.ncbi.nlm.nih.gov/37777343/)  
   <sub>The acne microbiome: from phyla to phylotypes — the most recent review.</sub>
33. Coenye T, Spittaels KJ, Achermann Y. The role of biofilm formation in the pathogenesis and antimicrobial susceptibility of Cutibacterium acnes. *Biofilm*. 2022;4:100063. [PMID 34950868](https://pubmed.ncbi.nlm.nih.gov/34950868/)  
   <sub>The role biofilm formation plays in pathogenesis and in antimicrobial susceptibility — the CAB compartment and the protection factor PROTB.</sub>
34. Ruffier d'Epenoux L, Fayoux E, Veziers J, et al. Biofilm of Cutibacterium acnes: a target of different active substances. *Int J Dermatol*. 2024;63:1541-1550. [PMID 38760974](https://pubmed.ncbi.nlm.nih.gov/38760974/)  
   <sub>A review of active substances that target the C. acnes biofilm.</sub>
35. Johnson T, Kang D, Barnard E, et al. Strain-Level Differences in Porphyrin Production and Regulation in Propionibacterium acnes Elucidate Disease Associations. *mSphere*. 2016;1. [PMID 27303708](https://pubmed.ncbi.nlm.nih.gov/27303708/)  
   <sub>Strain-level differences in porphyrin production and its regulation — the PORP state variable and its dependence on PHYLIA.</sub>
36. Nakatsuji T, Tang DC, Zhang L, et al. Propionibacterium acnes CAMP factor and host acid sphingomyelinase contribute to bacterial virulence: potential targets for inflammatory acne treatment. *PLoS One*. 2011;6:e14797. [PMID 21533261](https://pubmed.ncbi.nlm.nih.gov/21533261/)  
   <sub>CAMP factor and host acid sphingomyelinase — the keratinocyte cytotoxicity node.</sub>
37. Nakase K, Momose M, Yukawa T, et al. Development of skin sebum medium and inhibition of lipase activity in Cutibacterium acnes by oleic acid. *Access Microbiol*. 2022;4:acmi000397. [PMID 36415741](https://pubmed.ncbi.nlm.nih.gov/36415741/)  
   <sub>Development of a sebum medium and inhibition of C. acnes lipase activity — the lipase → free fatty acid (FFA) route.</sub>
38. Wang Y, Kao MS, Yu J, et al. A Precision Microbiome Approach Using Sucrose for Selective Augmentation of Staphylococcus epidermidis Fermentation against Propionibacterium acnes. *Int J Mol Sci*. 2016;17. [PMID 27834859](https://pubmed.ncbi.nlm.nih.gov/27834859/)  
   <sub>Inhibition of C. acnes using S. epidermidis fermentation — the antagonistic flora node.</sub>

## F. Innate immunity · the inflammasome · Th17 (innate and adaptive inflammation)

39. Kim J. Review of the innate immune response in acne vulgaris: activation of Toll-like receptor 2 in acne triggers inflammatory cytokine responses. *Dermatology*. 2005;211:193-8. [PMID 16205063](https://pubmed.ncbi.nlm.nih.gov/16205063/)  
   <sub>A review of the innate immune response in acne — TLR2 activation. The grounds for the state variable TLR.</sub>
40. Qin M, Pirouz A, Kim MH, et al. Propionibacterium acnes Induces IL-1β secretion via the NLRP3 inflammasome in human monocytes. *J Invest Dermatol*. 2014;134:381-388. [PMID 23884315](https://pubmed.ncbi.nlm.nih.gov/23884315/)  
   <sub>C. acnes makes human monocytes secrete IL-1β through the NLRP3 inflammasome — the IL1B production term and the NLRP3G second signal.</sub>
41. Nagy I, Pivarcsi A, Koreck A, et al. Distinct strains of Propionibacterium acnes induce selective human beta-defensin-2 and interleukin-8 expression in human keratinocytes through toll-like receptors. *J Invest Dermatol*. 2005;124:931-8. [PMID 15854033](https://pubmed.ncbi.nlm.nih.gov/15854033/)  
   <sub>hBD-2 and IL-8 expression is induced selectively according to the strain — phylotype-dependent inflammation.</sub>
42. Agak GW, Qin M, Nobe J, et al. Propionibacterium acnes Induces an IL-17 Response in Acne Vulgaris that Is Regulated by Vitamin A and Vitamin D. *J Invest Dermatol*. 2014;134:366-373. [PMID 23924903](https://pubmed.ncbi.nlm.nih.gov/23924903/)  
   <sub>C. acnes induces an IL-17 response, and vitamin A and D modulate it — the IL17 compartment.</sub>
43. Kelhälä HL, Palatsi R, Fyhrquist N, et al. IL-17/Th17 pathway is activated in acne lesions. *PLoS One*. 2014;9:e105238. [PMID 25153527](https://pubmed.ncbi.nlm.nih.gov/25153527/)  
   <sub>A human tissue study showing that the IL-17/Th17 pathway is genuinely activated in acne lesions.</sub>
44. Kang S, Cho S, Chung JH, et al. Inflammation and extracellular matrix degradation mediated by activated transcription factors nuclear factor-kappaB and activator protein-1 in inflammatory acne lesions in vivo. *Am J Pathol*. 2005;166:1691-9. [PMID 15920154](https://pubmed.ncbi.nlm.nih.gov/15920154/)  
   <sub>NF-κB and AP-1 activation and matrix metalloproteinases within acne lesions — the MMP state variable and the scarring route.</sub>

## G. Neuroendocrine · stress axis

45. Krause K, Schnitger A, Fimmel S, et al. Corticotropin-releasing hormone skin signaling is receptor-mediated and is predominant in the sebaceous glands. *Horm Metab Res*. 2007;39:166-70. [PMID 17326013](https://pubmed.ncbi.nlm.nih.gov/17326013/)  
   <sub>CRH signalling in skin is receptor-mediated and predominates in the sebaceous gland — the CRH → sebum production (GCRH) term.</sub>
46. Chiu A, Chon SY, Kimball AB. The response of skin disease to stress: changes in the severity of acne vulgaris as affected by examination stress. *Arch Dermatol*. 2003;139:897-900. [PMID 12873885](https://pubmed.ncbi.nlm.nih.gov/12873885/)  
   <sub>The change in acne severity with examination stress — the clinical grounds for the STRESS parameter.</sub>

## H. Scarring · post-inflammatory hyperpigmentation (PIH)

47. Connolly D, Vu HL, Mariwalla K, et al. Acne Scarring-Pathogenesis, Evaluation, and Treatment Options. *J Clin Aesthet Dermatol*. 2017;10:12-23. [PMID 29344322](https://pubmed.ncbi.nlm.nih.gov/29344322/)  
   <sub>A review of the pathogenesis, assessment and treatment of acne scarring — design of the SCAR state variable.</sub>
48. Tan J, Beissert S, Cook-Bolden F, et al. Evaluation of psychological well-being and social impact of atrophic acne scarring: A multinational, mixed-methods study. *JAAD Int*. 2022;6:43-50. [PMID 35005652](https://pubmed.ncbi.nlm.nih.gov/35005652/)  
   <sub>The psychological and social impact of atrophic acne scars.</sub>
49. Callender VD, St Surin-Lord S, Davis EC, et al. Postinflammatory hyperpigmentation: etiologic and therapeutic considerations. *Am J Clin Dermatol*. 2011;12:87-99. [PMID 21348540](https://pubmed.ncbi.nlm.nih.gov/21348540/)  
   <sub>The causes and treatment of post-inflammatory hyperpigmentation — the PIH state variable and the SKINPIG (Fitzpatrick) parameter.</sub>
50. Dréno B, Bissonnette R, Gagné-Henley A, et al. Prevention and Reduction of Atrophic Acne Scars with Adapalene 0.3%/Benzoyl Peroxide 2.5% Gel in Subjects with Moderate or Severe Facial Acne: Results of a 6-Month Randomized, Vehicle-Controlled Trial Using Intra-Individual Comparison. *Am J Clin Dermatol*. 2018;19:275-286. [PMID 29549588](https://pubmed.ncbi.nlm.nih.gov/29549588/)  
   <sub>Prevention and reduction of atrophic scarring by adapalene 0.3%/BPO — the clinical counterpart of the model's conclusion that 'the only way to reduce scarring is to reduce the time a nodule is present'.</sub>

## I. Topical therapy

51. Kolli SS, Pecone D, Pona A, et al. Topical Retinoids in Acne Vulgaris: A Systematic Review. *Am J Clin Dermatol*. 2019;20:345-365. [PMID 30674002](https://pubmed.ncbi.nlm.nih.gov/30674002/)  
   <sub>A systematic review of topical retinoids — the calibration target for the 12-week rate of lesion reduction.</sub>
52. Thiboutot DM, Weiss J, Bucko A, et al. Adapalene-benzoyl peroxide, a fixed-dose combination for the treatment of acne vulgaris: results of a multicenter, randomized double-blind, controlled study. *J Am Acad Dermatol*. 2007;57:791-9. [PMID 17655969](https://pubmed.ncbi.nlm.nih.gov/17655969/)  
   <sub>The adapalene–BPO fixed-dose combination — the grounds for scenario 5.</sub>
53. Del Rosso JQ, Lain E, York JP, et al. Trifarotene 0.005% Cream in the Treatment of Facial and Truncal Acne Vulgaris in Patients with Skin of Color: a Case Series. *Dermatol Ther (Heidelb)*. 2022;12:2189-2200. [PMID 35994159](https://pubmed.ncbi.nlm.nih.gov/35994159/)  
   <sub>The trifarotene 0.005% cream trials in facial and truncal acne — sets the relative potency in ACN_RETPOT.</sub>
54. Tanghetti E, Dhawan S, Green L, et al. Randomized comparison of the safety and efficacy of tazarotene 0.1% cream and adapalene 0.3% gel in the treatment of patients with at least moderate facial acne vulgaris. *J Drugs Dermatol*. 2010;9:549-58. [PMID 20480800](https://pubmed.ncbi.nlm.nih.gov/20480800/)  
   <sub>A randomised comparison of tazarotene 0.1% cream and adapalene.</sub>
55. Fulton JE Jr, Farzad-Bakshandeh A, Bradley S. Studies on the mechanism of action to topical benzoyl peroxide and vitamin A acid in acne vulgaris. *J Cutan Pathol*. 1974;1:191-200. [PMID 4283462](https://pubmed.ncbi.nlm.nih.gov/4283462/)  
   <sub>The classic study of the mechanism of action of topical benzoyl peroxide and vitamin A acid.</sub>
56. Leyden JJ, Preston N, Osborn C, et al. In-vivo Effectiveness of Adapalene 0.1%/Benzoyl Peroxide 2.5% Gel on Antibiotic-sensitive and Resistant Propionibacterium acnes. *J Clin Aesthet Dermatol*. 2011;4:22-6. [PMID 21607190](https://pubmed.ncbi.nlm.nih.gov/21607190/)  
   <sub>The in vivo effect of adapalene 0.1%/BPO 2.5% gel against antibiotic-resistant P. acnes — the direct grounds for the KBPORES term, by which BPO clears the resistant fraction.</sub>
57. Webster G, Thiboutot DM, Chen DM, et al. Impact of a fixed combination of clindamycin phosphate 1.2%-benzoyl peroxide 2.5% aqueous gel on health-related quality of life in moderate to severe acne vulgaris. *Cutis*. 2010;86:263-7. [PMID 21214129](https://pubmed.ncbi.nlm.nih.gov/21214129/)  
   <sub>The effect of the clindamycin 1.2%–BPO fixed combination — scenario 7.</sub>
58. Gold M, Lain T, Harper JC, et al. Efficacy and Safety of Clindamycin Phosphate 1.2%/Adapalene 0.15%/Benzoyl Peroxide 3.1% Gel: Post Hoc Analysis by Baseline Disease Severity. *Dermatol Ther (Heidelb)*. 2025;15:1867-1882. [PMID 40377868](https://pubmed.ncbi.nlm.nih.gov/40377868/)  
   <sub>Efficacy and safety of the clindamycin 1.2%/adapalene 0.15%/BPO 3.1% triple fixed combination.</sub>
59. Draelos ZD, Rodriguez DA, Kempers SE, et al. Treatment Response With Once-Daily Topical Dapsone Gel, 7.5% for Acne Vulgaris: Subgroup Analysis of Pooled Data from Two Randomized, Double-Blind Stu. *J Drugs Dermatol*. 2017;16:591-598. [PMID 28686777](https://pubmed.ncbi.nlm.nih.gov/28686777/)  
   <sub>Dapsone gel 7.5% once daily — the neutrophil (MPO) inhibition term (EDAP).</sub>
60. Gowda CM, Wairkar S. Azelaic acid-based lyotropic liquid crystals gel for acne vulgaris: Formulation optimization, antimicrobial activity and dermatopharmacokinetic study. *Int J Pharm*. 2024;667:124879. [PMID 39490554](https://pubmed.ncbi.nlm.nih.gov/39490554/)  
   <sub>A study of azelaic acid formulations — the multiple actions of normalising keratinisation, antimicrobial activity and tyrosinase inhibition.</sub>
61. Makrantonaki E. Topical minocycline foam: A new option for acne treatment?. *J Eur Acad Dermatol Venereol*. 2025;39:885-886. [PMID 40277219](https://pubmed.ncbi.nlm.nih.gov/40277219/)  
   <sub>Topical minocycline foam — a topical antibiotic option that minimises systemic exposure.</sub>
62. Stein Gold L, Kwong P, Draelos Z, et al. Impact of Topical Vehicles and Cutaneous Delivery Technologies on Patient Adherence and Treatment Outcomes in Acne and Rosacea. *J Clin Aesthet Dermatol*. 2023;16:26-34. [PMID 37288283](https://pubmed.ncbi.nlm.nih.gov/37288283/)  
   <sub>The influence of vehicle and delivery technology on adherence — the ADHERE parameter and the irritation-adherence feedback.</sub>

## J. Systemic antibiotics · antibiotic resistance and stewardship

63. Garner SE, Eady A, Bennett C, et al. Minocycline for acne vulgaris: efficacy and safety. *Cochrane Database Syst Rev*. 2012;2012:CD002086. [PMID 22895927](https://pubmed.ncbi.nlm.nih.gov/22895927/)  
   <sub>The efficacy and safety of minocycline (Cochrane) — no evidence of superiority within the class.</sub>
64. Zhanel G, Critchley I, Lin LY, et al. Microbiological Profile of Sarecycline, a Novel Targeted Spectrum Tetracycline for the Treatment of Acne Vulgaris. *Antimicrob Agents Chemother*. 2019;63. [PMID 30397052](https://pubmed.ncbi.nlm.nih.gov/30397052/)  
   <sub>The microbiological profile of sarecycline — the narrow-spectrum (NARROW) parameter.</sub>
65. Zhang J, He L, Chen X, et al. Efficacy and Safety of Sarecycline in Chinese Patients with Moderate-to-Severe Acne Vulgaris: Randomized Phase 3 Clinical Trial with Open-Label Follow-Up. *Dermatol Ther (Heidelb)*. 2025;15:3285-3300. [PMID 40877730](https://pubmed.ncbi.nlm.nih.gov/40877730/)  
   <sub>The efficacy and safety of sarecycline in moderate–severe acne.</sub>
66. Toossi P, Farshchian M, Malekzad F, et al. Subantimicrobial-dose doxycycline in the treatment of moderate facial acne. *J Drugs Dermatol*. 2008;7:1149-52. [PMID 19137768](https://pubmed.ncbi.nlm.nih.gov/19137768/)  
   <sub>Sub-antimicrobial dose doxycycline in the treatment of moderate facial acne — the key evidence that there is an effect even without the antimicrobial arm.</sub>
67. Moore A, Ling M, Bucko A, et al. Efficacy and Safety of Subantimicrobial Dose, Modified-Release Doxycycline 40 mg Versus Doxycycline 100 mg Versus Placebo for the treatment of Inflammatory Lesions in Moderate and Severe Acne: A Randomized, Double-Blinded, Controlled Study. *J Drugs Dermatol*. 2015;14:581-6. [PMID 26091383](https://pubmed.ncbi.nlm.nih.gov/26091383/)  
   <sub>The efficacy and safety of sustained-release sub-antimicrobial dose doxycycline 40 mg — the grounds for a design that separates the two arms as EC5AI ≪ EC5TET.</sub>
68. Ross JI, Snelling AM, Eady EA, et al. Phenotypic and genotypic characterization of antibiotic-resistant Propionibacterium acnes isolated from acne patients attending dermatology clinics in Europe, the U.S.A., Japan and Australia. *Br J Dermatol*. 2001;144:339-46. [PMID 11251569](https://pubmed.ncbi.nlm.nih.gov/11251569/)  
   <sub>Phenotypic and genotypic characteristics of antibiotic-resistant P. acnes recovered from acne patients.</sub>
69. Ross JI, Eady EA, Cove JH, et al. 16S rRNA mutation associated with tetracycline resistance in a gram-positive bacterium. *Antimicrob Agents Chemother*. 1998;42:1702-5. [PMID 9661007](https://pubmed.ncbi.nlm.nih.gov/9661007/)  
   <sub>The 16S rRNA mutation associated with tetracycline resistance — the molecular substance of RESF.</sub>
70. Zhu C, Wei B, Li Y, et al. Antibiotic resistance rates in Cutibacterium acnes isolated from patients with acne vulgaris: a systematic review and meta-analysis. *Front Microbiol*. 2025;16:1565111. [PMID 40535003](https://pubmed.ncbi.nlm.nih.gov/40535003/)  
   <sub>The most recent survey of antibiotic resistance rates in C. acnes.</sub>
71. Barbieri JS, Spaccarelli N, Margolis DJ, et al. Approaches to limit systemic antibiotic use in acne: Systemic alternatives, emerging topical therapies, dietary modification, and laser and light-based treatments. *J Am Acad Dermatol*. 2019;80:538-549. [PMID 30296534](https://pubmed.ncbi.nlm.nih.gov/30296534/)  
   <sub>Approaches to limiting systemic antibiotic use in acne.</sub>
72. Rosenberg AL, Shah M, Del Rosso JQ, et al. Optimal Use Recommendations and Stewardship Principles with Oral Antibiotics in Acne Vulgaris Management: An Expert Consensus Panel. *J Clin Aesthet Dermatol*. 2025;18:21-29. [PMID 41640785](https://pubmed.ncbi.nlm.nih.gov/41640785/)  
   <sub>Recommendations for the optimal use of oral antibiotics and the principles of stewardship — the clinical counterpart of the scenario 6 vs 7 comparison.</sub>

## K. Hormonal therapy

73. Arowojolu AO, Gallo MF, Lopez LM, et al. Combined oral contraceptive pills for treatment of acne. *Cochrane Database Syst Rev*. 2012;2012:CD004425. [PMID 22786490](https://pubmed.ncbi.nlm.nih.gov/22786490/)  
   <sub>The effect of combined oral contraceptives in the treatment of acne (Cochrane) — the effect size for scenario 12.</sub>
74. Koo EB, Petersen TD, Kimball AB. Meta-analysis comparing efficacy of antibiotics versus oral contraceptives in acne vulgaris. *J Am Acad Dermatol*. 2014;71:450-9. [PMID 24880665](https://pubmed.ncbi.nlm.nih.gov/24880665/)  
   <sub>A meta-analysis comparing the efficacy of antibiotics against oral contraceptives — equivalent at 6 months.</sub>
75. Santer M, Lawrence M, Pyne S, et al. Clinical and cost-effectiveness of spironolactone in treating persistent facial acne in women: SAFA double-blinded RCT. *Health Technol Assess*. 2024;28:1-86. [PMID 39268864](https://pubmed.ncbi.nlm.nih.gov/39268864/)  
   <sub>The SAFA randomised double-blind trial — the clinical and cost effectiveness of spironolactone in adult female acne.</sub>
76. Kow CS, Ramachandram DS, Hasan SS, et al. Spironolactone for the Treatment of Moderate to Severe Acne in Adult Women: A Systematic Review and Meta-Analysis of Randomised Controlled Trials. *Australas J Dermatol*. 2025;66:165-168. [PMID 39912292](https://pubmed.ncbi.nlm.nih.gov/39912292/)  
   <sub>Spironolactone in moderate–severe acne in adult women.</sub>
77. Basu P, Elman SA, Abudu B, et al. High-dose spironolactone for acne in patients with polycystic ovarian syndrome: A single-institution retrospective study. *J Am Acad Dermatol*. 2021;85:740-741. [PMID 31400460](https://pubmed.ncbi.nlm.nih.gov/31400460/)  
   <sub>High-dose spironolactone in patients with PCOS — the pcos phenotype scenario.</sub>
78. Patiyasikunt M, Chancheewa B, Asawanonda P, et al. Efficacy and tolerability of low-dose spironolactone and topical benzoyl peroxide in adult female acne: A randomized, double-blind, placebo-controlled trial. *J Dermatol*. 2020;47:1411-1416. [PMID 32857471](https://pubmed.ncbi.nlm.nih.gov/32857471/)  
   <sub>The efficacy and tolerability of low-dose spironolactone combined with topical BPO.</sub>
79. Plovanich M, Weng QY, Mostaghimi A. Low Usefulness of Potassium Monitoring Among Healthy Young Women Taking Spironolactone for Acne. *JAMA Dermatol*. 2015;151:941-4. [PMID 25796182](https://pubmed.ncbi.nlm.nih.gov/25796182/)  
   <sub>The low utility of potassium monitoring in healthy young women — the interpretation of the KSER tracking state variable.</sub>
80. Rosette C, Agan FJ, Mazzetti A, et al. Cortexolone 17α-propionate (Clascoterone) Is a Novel Androgen Receptor Antagonist that Inhibits Production of Lipids and Inflammatory Cytokines from Sebocytes In Vitro. *J Drugs Dermatol*. 2019;18:412-418. [PMID 31141847](https://pubmed.ncbi.nlm.nih.gov/31141847/)  
   <sub>Clascoterone (cortexolone 17α-propionate), as a sebocyte AR antagonist, suppresses lipid and inflammatory cytokine production — the KICLAS competitive antagonism term.</sub>
81. Alkhodaidi ST, Al Hawsawi KA, Alkhudaidi IT, et al. Efficacy and safety of topical clascoterone cream for treatment of acne vulgaris: A systematic review and meta-analysis of randomized placebo-controlled trials. *Dermatol Ther*. 2021;34:e14609. [PMID 33258536](https://pubmed.ncbi.nlm.nih.gov/33258536/)  
   <sub>The efficacy and safety of topical clascoterone cream — scenario 11.</sub>
82. Damoulaki E, Sioutis D, Sarli V, et al. Polycystic Ovary Syndrome-Associated Acne: The Interplay of Hyperandrogenism, Insulin Resistance, and Therapeutic Strategies. *Cureus*. 2025;17:e98103. [PMID 41473651](https://pubmed.ncbi.nlm.nih.gov/41473651/)  
   <sub>The interaction between PCOS-associated acne and hyperandrogenaemia.</sub>

## L. Isotretinoin — PK · mechanism of action · cumulative dose · safety

83. Colburn WA, Gibson DM, Wiens RE, et al. Food increases the bioavailability of isotretinoin. *J Clin Pharmacol*. 1983;23:534-9. [PMID 6582073](https://pubmed.ncbi.nlm.nih.gov/6582073/)  
   <sub>Food raises the bioavailability of isotretinoin — the original grounds for the FOOD/FOODEF parameters.</sub>
84. Webster GF, Leyden JJ, Gross JA. Results of a Phase III, double-blind, randomized, parallel-group, non-inferiority study evaluating the safety and efficacy of isotretinoin-Lidose in patients with severe recalcitrant nodular acne. *J Drugs Dermatol*. 2014;13:665-70. [PMID 24918555](https://pubmed.ncbi.nlm.nih.gov/24918555/)  
   <sub>The phase 3 non-inferiority trial of Lidose isotretinoin — the LIDOSE flag (loss of food dependence).</sub>
85. Almond-Roesler B, Blume-Peytavi U, Bisson S, et al. Monitoring of isotretinoin therapy by measuring the plasma levels of isotretinoin and 4-oxo-isotretinoin. A useful tool for management of severe acne. *Dermatology*. 1998;196:176-81. [PMID 9557257](https://pubmed.ncbi.nlm.nih.gov/9557257/)  
   <sub>Monitoring of plasma concentrations of isotretinoin and 4-oxo-isotretinoin — the metabolite compartment and POTOXO.</sub>
86. Landthaler M, Kummermehr J, Wagner A, et al. Inhibitory effects of 13-cis-retinoic acid on human sebaceous glands. *Arch Dermatol Res*. 1980;269:297-309. [PMID 6453562](https://pubmed.ncbi.nlm.nih.gov/6453562/)  
   <sub>The suppressive effect of 13-cis-retinoic acid on the human sebaceous gland — the target of a 90% reduction in SER.</sub>
87. Nelson AM, Gilliland KL, Cong Z, et al. 13-cis Retinoic acid induces apoptosis and cell cycle arrest in human SEB-1 sebocytes. *J Invest Dermatol*. 2006;126:2178-89. [PMID 16575387](https://pubmed.ncbi.nlm.nih.gov/16575387/)  
   <sub>13-cis retinoic acid induces apoptosis and cell cycle arrest in human SEB-1 sebocytes — the KAPO augmentation term (EISOAP).</sub>
88. Nelson AM, Zhao W, Gilliland KL, et al. Neutrophil gelatinase-associated lipocalin mediates 13-cis retinoic acid-induced apoptosis of human sebaceous gland cells. *J Clin Invest*. 2008;118:1468-78. [PMID 18317594](https://pubmed.ncbi.nlm.nih.gov/18317594/)  
   <sub>NGAL (lipocalin-2) mediates 13-cis retinoic acid-induced sebocyte apoptosis — the molecular route of cluster 15 of the map.</sub>
89. Layton AM, Knaggs H, Taylor J, et al. Isotretinoin for acne vulgaris--10 years later: a safe and successful treatment. *Br J Dermatol*. 1993;129:292-6. [PMID 8286227](https://pubmed.ncbi.nlm.nih.gov/8286227/)  
   <sub>10 years of experience with isotretinoin — the clinical determinants of relapse.</sub>
90. Blasiak RC, Stamey CR, Burkhart CN, et al. High-dose isotretinoin treatment and the rate of retrial, relapse, and adverse effects in patients with acne vulgaris. *JAMA Dermatol*. 2013;149:1392-8. [PMID 24173086](https://pubmed.ncbi.nlm.nih.gov/24173086/)  
   <sub>High-dose isotretinoin and the rates of repeat course, relapse and adverse events.</sub>
91. Borghi A, Mantovani L, Minghetti S, et al. Low-cumulative dose isotretinoin treatment in mild-to-moderate acne: efficacy in achieving stable remission. *J Eur Acad Dermatol Venereol*. 2011;25:1094-8. [PMID 21198947](https://pubmed.ncbi.nlm.nih.gov/21198947/)  
   <sub>The ability of low cumulative dose isotretinoin to achieve stable remission — the comparator for scenario 15.</sub>
92. Rademaker M. Making sense of the effects of the cumulative dose of isotretinoin in acne vulgaris. *Int J Dermatol*. 2016;55:518-23. [PMID 26471145](https://pubmed.ncbi.nlm.nih.gov/26471145/)  
   <sub>Interpreting the cumulative dose effect of isotretinoin in acne — the grounds for the CUMISO → DURAB Hill function (CD50 ≈ 85 mg/kg).</sub>
93. Zane LT, Leyden WA, Marqueling AL, et al. A population-based analysis of laboratory abnormalities during isotretinoin therapy for acne vulgaris. *Arch Dermatol*. 2006;142:1016-22. [PMID 16924051](https://pubmed.ncbi.nlm.nih.gov/16924051/)  
   <sub>A population-based analysis of laboratory abnormalities during isotretinoin treatment — calibration of the TG and ALT state variables.</sub>
94. Xia E, Han J, Faletsky A, et al. Isotretinoin Laboratory Monitoring in Acne Treatment: A Delphi Consensus Study. *JAMA Dermatol*. 2022;158:942-948. [PMID 35704293](https://pubmed.ncbi.nlm.nih.gov/35704293/)  
   <sub>A Delphi consensus on laboratory monitoring for isotretinoin.</sub>
95. Huang YC, Cheng YC. Isotretinoin treatment for acne and risk of depression: A systematic review and meta-analysis. *J Am Acad Dermatol*. 2017;76:1068-1076.e9. [PMID 28291553](https://pubmed.ncbi.nlm.nih.gov/28291553/)  
   <sub>A systematic review and meta-analysis of isotretinoin and the risk of depression.</sub>
96. Coberly S, Lammer E, Alashari M. Retinoic acid embryopathy: case report and review of literature. *Pediatr Pathol Lab Med*. 1996;16:823-36. [PMID 9025880](https://pubmed.ncbi.nlm.nih.gov/9025880/)  
   <sub>Retinoic acid embryopathy — the grounds for teratogenicity and for the iPLEDGE requirements.</sub>
97. Fraunfelder FW, Fraunfelder FT, Corbett JJ. Isotretinoin-associated intracranial hypertension. *Ophthalmology*. 2004;111:1248-50. [PMID 15177980](https://pubmed.ncbi.nlm.nih.gov/15177980/)  
   <sub>Isotretinoin-associated raised intracranial pressure — the grounds for the contraindication to concomitant tetracycline.</sub>
98. Scaramuzzino L, Coronella L, Lauletta G, et al. Recalcitrant isotretinoin-induced acne fulminans successfully treated with oral dapsone. *Ital J Dermatol Venerol*. 2025;160:381-382. [PMID 40292612](https://pubmed.ncbi.nlm.nih.gov/40292612/)  
   <sub>Isotretinoin-induced acne fulminans — the initial flare in severe nodular disease and the need for a preceding steroid.</sub>

## M. Clinical practice guidelines

99. Reynolds RV, Yeung H, Cheng CE, et al. Guidelines of care for the management of acne vulgaris. *J Am Acad Dermatol*. 2024;90:1006.e1-1006.e30. [PMID 38300170](https://pubmed.ncbi.nlm.nih.gov/38300170/)  
   <sub>The American Academy of Dermatology (AAD) 2024 acne guideline — the basis on which regimens are composed in the scenario library.</sub>
100. Zaenglein AL, Pathy AL, Schlosser BJ, et al. Guidelines of care for the management of acne vulgaris. *J Am Acad Dermatol*. 2016;74:945-73.e33. [PMID 26897386](https://pubmed.ncbi.nlm.nih.gov/26897386/)  
   <sub>The AAD 2016 acne guideline.</sub>
101. Nast A, Al Wattar BH, Beylot Barry M, et al. Update of the EuroGuiDerm evidence-based guideline for the treatment of acne-Short version. *J Eur Acad Dermatol Venereol*. 2026;40:1162-1172. [PMID 41847993](https://pubmed.ncbi.nlm.nih.gov/41847993/)  
   <sub>The EuroGuiDerm evidence-based acne treatment guideline update.</sub>
102. Nast A, Rosumeck S, Sammain A, et al. Methods report on the development of the European S3 guidelines for the treatment of acne. *J Eur Acad Dermatol Venereol*. 2012;26 Suppl 1:e1-41. [PMID 22356612](https://pubmed.ncbi.nlm.nih.gov/22356612/)  
   <sub>A report of the methodology used to develop the European S3 guideline.</sub>
103. Tan J, Alexis A, Baldwin H, et al. The Personalised Acne Care Pathway-Recommendations to guide longitudinal management from the Personalising Acne: Consensus of Experts. *JAAD Int*. 2021;5:101-111. [PMID 34816135](https://pubmed.ncbi.nlm.nih.gov/34816135/)  
   <sub>Recommendations for a personalised acne management pathway — treatment choice by phenotype.</sub>
104. Thiboutot DM, Shalita AR, Yamauchi PS, et al. Adapalene gel, 0.1%, as maintenance therapy for acne vulgaris: a randomized, controlled, investigator-blind follow-up of a recent combination study. *Arch Dermatol*. 2006;142:597-602. [PMID 16702497](https://pubmed.ncbi.nlm.nih.gov/16702497/)  
   <sub>A randomised trial of adapalene 0.1% gel as maintenance therapy — the direct grounds for scenario 17 (maintenance after induction).</sub>

## N. Variant phenotypes · differential diagnosis

105. Gorji M, Joseph J, Pavlakis N, et al. Prevention and management of acneiform rash associated with EGFR inhibitor therapy: A systematic review and meta-analysis. *Asia Pac J Clin Oncol*. 2022;18:526-539. [PMID 35352492](https://pubmed.ncbi.nlm.nih.gov/35352492/)  
   <sub>Prevention and management of the acneiform eruption associated with EGFR inhibitors.</sub>
106. Wang M, Zhang P, Shen M. A de novo heterozygous PSTPIP1 variant associated with PAPA syndrome: a Chinese case report and literature review. *Front Genet*. 2026;17:1825761. [PMID 42358432](https://pubmed.ncbi.nlm.nih.gov/42358432/)  
   <sub>PSTPIP1 mutation and PAPA syndrome — the autoinflammatory overlap phenotype.</sub>

## O. QSP methodology · implementation

107. Kapitanov GI, Earp JC, Gadkar K, et al. Bridging the Gap: Integrating Quantitative Systems Pharmacology and Pharmacometrics in Drug Development. *Clin Pharmacol Ther*. 2026;119:830-833. [PMID 41472478](https://pubmed.ncbi.nlm.nih.gov/41472478/)  
   <sub>The integration of QSP and pharmacometrics.</sub>
108. Cheng Y, Straube R, Alnaif AE, et al. Virtual Populations for Quantitative Systems Pharmacology Models. *Methods Mol Biol*. 2022;2486:129-179. [PMID 35437722](https://pubmed.ncbi.nlm.nih.gov/35437722/)  
   <sub>Building a virtual population for a QSP model — the methodological background to this model's P_success (logistic transform) approach.</sub>
109. Coto-Segura P, Segú-Vergés C, Martorell A, et al. A quantitative systems pharmacology model for certolizumab pegol treatment in moderate-to-severe psoriasis. *Front Immunol*. 2023;14:1212981. [PMID 37809085](https://pubmed.ncbi.nlm.nih.gov/37809085/)  
   <sub>An example of a QSP model in a skin disease (psoriasis) — a precedent for dermatological QSP.</sub>
110. Lu T, Poon V, Brooks L, et al. gPKPDviz: A flexible R shiny tool for pharmacokinetic/pharmacodynamic simulations using mrgsolve. *CPT Pharmacometrics Syst Pharmacol*. 2024;13:341-358. [PMID 38082557](https://pubmed.ncbi.nlm.nih.gov/38082557/)  
   <sub>gPKPDviz — an mrgsolve-based Shiny tool for PK/PD simulation.</sub>

---

## Where the literature enters the model structure (mapping)

| Model element | Supporting literature (section) |
|---|---|
| The four pathogenic pillars and the microcomedone reservoir (`MC`) | A · B |
| Androgen → `ARS` → sebaceous gland mass `SGM` · sebum rate `SER` | C |
| IGF-1 / insulin / mTORC1 · FoxO1 → `LIP` | D |
| `SER` rise → linoleic acid dilution `LA` → hyperkeratinisation `KER` | B · C |
| Squalene peroxide `SQOX` (comedogenic + NLRP3 second signal) | C · F |
| `CAP`/`CAB` logistic growth, niche = sebum, biofilm protection | E |
| Resistant fraction `RESF` — antibiotic selection pressure ↑ / BPO clearance ↓ | J |
| TLR2 → NLRP3 → IL-1β → IL-8 → neutrophil → MMP | F |
| Th17/IL-17 amplification loop | F |
| CRH·stress → sebum production | G |
| The lesion transition chain and the accumulation of scarring·PIH | B · H |
| Topical retinoid/BPO/clindamycin/dapsone/azelaic acid PD | I |
| Separation of the antimicrobial arm and the non-antimicrobial anti-inflammatory arm of tetracycline | J |
| COC(SHBG·LH) · spironolactone · clascoterone | K |
| Isotretinoin PK(food effect·4-oxo), sebaceous gland apoptosis, cumulative dose–relapse | L |
| Safety(TG·ALT·K⁺·mucocutaneous·teratogenicity·IIH) | K · L |
| Scenario design and endpoint definition | I · J · K · M |
| QSP methodology·mrgsolve implementation | O |

---

**110 in total.** Every link points to PubMed.
This is model documentation for educational and research purposes; it does not replace clinical practice guidelines.

