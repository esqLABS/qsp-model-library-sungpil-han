# Urea Cycle Disorders (UCD) — QSP model references
### References for the UCD / OTC-deficiency Quantitative Systems Pharmacology model

This document collects the evidence for what went into `ucd_qsp_model.dot` (the mechanistic map), `ucd_mrgsolve_model.R` (the ODE model),
and `ucd_shiny_app.R` (the dashboard): **the basis for its structure, its parameters, and its clinical anchors**.
The head of each section records **which part of the model** that group of references fixed.

---

## 1. Reviews · guidelines

> Fixes the overall structure of the model, the diagnostic algorithm (citrulline + orotic acid discrimination), the ammonia thresholds
> (100 / 200 / 360 µmol/L), and the order of emergency management.

1. Häberle J, et al. **Suggested guidelines for the diagnosis and management of urea cycle disorders: first revision.** *J Inherit Metab Dis* 2019;42:1192-1230. — the revised European guideline. Primary basis for the ammonia thresholds, the indications for dialysis, and scavenger dosing. <https://pubmed.ncbi.nlm.nih.gov/30982989/>
2. Häberle J, et al. **Suggested guidelines for the diagnosis and management of urea cycle disorders.** *Orphanet J Rare Dis* 2012;7:32. <https://pubmed.ncbi.nlm.nih.gov/22642880/>
3. Summar ML, et al. **The incidence of urea cycle disorders.** *Mol Genet Metab* 2013;110:179-180. <https://pubmed.ncbi.nlm.nih.gov/23972786/>
4. Brusilow SW, Maestri NE. **Urea cycle disorders: diagnosis, pathophysiology, and therapy.** *Adv Pediatr* 1996;43:127-170. — the origin of the alternative nitrogen excretion pathway concept. <https://pubmed.ncbi.nlm.nih.gov/8794176/>
5. Matsumoto S, et al. **Urea cycle disorders — update.** *J Hum Genet* 2019;64:833-847. <https://pubmed.ncbi.nlm.nih.gov/31110235/>
6. Ah Mew N, et al. **Urea Cycle Disorders Overview.** *GeneReviews* (updated 2017). <https://pubmed.ncbi.nlm.nih.gov/20301396/>
7. Batshaw ML, Tuchman M, Summar M, Seminara J; Members of the UCDC. **A longitudinal study of urea cycle disorders.** *Mol Genet Metab* 2014;113:127-130. — the Urea Cycle Disorders Consortium longitudinal cohort. <https://pubmed.ncbi.nlm.nih.gov/25135652/>
8. Stepien KM, et al. **Challenges in the management of adult patients with urea cycle disorders.** *J Clin Med* 2022;11:2029. <https://pubmed.ncbi.nlm.nih.gov/35407637/>

## 2. Urea cycle enzymology · NAG-CPS1 regulation

> Fixes `VMAX_UC` (2400 µmol N/kg/h), `KMNH4` (120 µmol/L), the NAG activation term `f_NAG`,
> the stimulation of NAGS by arginine (`GARG`), and the activation of CPS1 by carglumic acid (`GNCG`).

9. Meijer AJ, Lamers WH, Chamuleau RA. **Nitrogen metabolism and ornithine cycle function.** *Physiol Rev* 1990;70:701-748. — the standard review of periportal/perivenous compartmentalisation in the liver and of urea cycle regulation. <https://pubmed.ncbi.nlm.nih.gov/2194222/>
10. Morris SM Jr. **Regulation of enzymes of the urea cycle and arginine metabolism.** *Annu Rev Nutr* 2002;22:87-105. <https://pubmed.ncbi.nlm.nih.gov/12055339/>
11. Caldovic L, Tuchman M. **N-acetylglutamate and its changing role through evolution.** *Biochem J* 2003;372:279-290. <https://pubmed.ncbi.nlm.nih.gov/12633501/>
12. Shi D, Allewell NM, Tuchman M. **The N-acetylglutamate synthase family: structures, function and mechanisms.** *Int J Mol Sci* 2015;16:13004-13022. <https://pubmed.ncbi.nlm.nih.gov/26068232/>
13. Nissim I, et al. **Effects of a glucokinase activator on hepatic intermediary metabolism: study with 13C-isotopomer-based metabolomics.** *Biochem J* 2012;444:537-551. — quantification of hepatic ureagenesis flux. <https://pubmed.ncbi.nlm.nih.gov/22448977/>
14. Rubio V, Grisolía S. **Human carbamoylphosphate synthetase I.** *Enzyme* 1981;26:233-239. <https://pubmed.ncbi.nlm.nih.gov/7318135/>
15. Yudkoff M, et al. **In vivo measurement of ureagenesis with stable isotopes.** *J Inherit Metab Dis* 1998;21(Suppl 1):21-29. — the method for measuring ureagenesis flux with 15N/13C; the primary pharmacodynamic index in the gene therapy trials. <https://pubmed.ncbi.nlm.nih.gov/9686341/>

## 3. OTC deficiency — genetics, X-inactivation, the phenotype spectrum

> The meaning of the `OTCACT` parameter, hepatic mosaicism in heterozygous women (the `Xinact` node),
> and the basis for the neonatal versus late-onset distinction (`Onset_class`).

16. Tuchman M, et al. **Cross-sectional multicenter study of patients with urea cycle disorders in the United States.** *Mol Genet Metab* 2008;94:397-402. <https://pubmed.ncbi.nlm.nih.gov/18562231/>
17. Yamaguchi S, et al. **Mutations and polymorphisms in the human ornithine transcarbamylase (OTC) gene.** *Hum Mutat* 2006;27:626-632. <https://pubmed.ncbi.nlm.nih.gov/16786505/>
18. Maestri NE, Brusilow SW, Clissold DB, Bassett SS. **Long-term treatment of girls with ornithine transcarbamylase deficiency.** *N Engl J Med* 1996;335:855-859. <https://pubmed.ncbi.nlm.nih.gov/8778603/>
19. Batshaw ML, Msall M, Beaudet AL, Trojak J. **Risk of serious illness in heterozygotes for ornithine transcarbamylase deficiency.** *J Pediatr* 1986;108:236-241. <https://pubmed.ncbi.nlm.nih.gov/3944709/>
20. Gyato K, Wray J, Huang ZJ, Yudkoff M, Batshaw ML. **Metabolic and neuropsychological phenotype in women heterozygous for ornithine transcarbamylase deficiency.** *Ann Neurol* 2004;55:80-86. — executive function deficits even in 'asymptomatic' carriers. The ADHD-like phenotype node of the model. <https://pubmed.ncbi.nlm.nih.gov/14705115/>
21. Lichter-Konecki U, Caldovic L, Morizono H, Simpson K. **Ornithine Transcarbamylase Deficiency.** *GeneReviews* (updated 2022). <https://pubmed.ncbi.nlm.nih.gov/24006547/>

## 4. Ammonia neurotoxicity — the astrocytic glutamine osmolyte hypothesis

> The basis for the most important CNS module in the model (`GLNB` · `MINSB` · `BRWAT` · `ICP`).
> The `FCOMP` interaction — that myo-inositol is already depleted in chronic hyperammonaemia, so that there is no
> buffer left when an acute rise comes — comes from here.

22. Brusilow SW, Koehler RC, Traystman RJ, Cooper AJL. **Astrocyte glutamine synthetase: importance in hyperammonemic syndromes and potential target for therapy.** *Neurotherapeutics* 2010;7:452-470. — the definitive statement of the osmolyte hypothesis. <https://pubmed.ncbi.nlm.nih.gov/20880508/>
23. Butterworth RF. **Pathophysiology of hepatic encephalopathy: a new look at ammonia.** *Metab Brain Dis* 2002;17:221-227. <https://pubmed.ncbi.nlm.nih.gov/12602499/>
24. Norenberg MD, Rao KV, Jayakumar AR. **Mechanisms of ammonia-induced astrocyte swelling.** *Metab Brain Dis* 2005;20:303-318. — the mitochondrial permeability transition (mPTP) and the ROS/RNS pathway. <https://pubmed.ncbi.nlm.nih.gov/16382342/>
25. Häussinger D, Laubenberger J, vom Dahl S, et al. **Proton magnetic resonance spectroscopy studies on human brain myo-inositol in hypo-osmolarity and hepatic encephalopathy.** *Gastroenterology* 1994;107:1475-1480. — direct evidence of the reciprocal exchange glutamine↑ / myo-inositol↓. <https://pubmed.ncbi.nlm.nih.gov/7926510/>
26. Connelly A, Cross JH, Gadian DG, Hunter JV, Kirkham FJ, Leonard JV. **Magnetic resonance spectroscopy shows increased brain glutamine in ornithine carbamoyl transferase deficiency.** *Pediatr Res* 1993;33:77-81. <https://pubmed.ncbi.nlm.nih.gov/8433886/>
27. Gropman AL, et al. **1H MRS identifies symptomatic and asymptomatic subjects with partial ornithine transcarbamylase deficiency.** *Mol Genet Metab* 2008;95:21-30. <https://pubmed.ncbi.nlm.nih.gov/18662894/>
28. Gropman AL, Summar M, Leonard JV. **Neurological implications of urea cycle disorders.** *J Inherit Metab Dis* 2007;30:865-879. <https://pubmed.ncbi.nlm.nih.gov/18038189/>
29. Cooper AJL, Plum F. **Biochemistry and physiology of brain ammonia.** *Physiol Rev* 1987;67:440-519. — NH3/NH4+ partitioning, ion trapping, and dependence on brain pH. The basis for the `PKANH3` and `PHBR` parameters. <https://pubmed.ncbi.nlm.nih.gov/2882529/>
30. Bachmann C. **Mechanisms of hyperammonemia.** *Clin Chem Lab Med* 2002;40:653-662. <https://pubmed.ncbi.nlm.nih.gov/12241009/>
31. Cagnon L, Braissant O. **Hyperammonemia-induced toxicity for the developing central nervous system.** *Brain Res Rev* 2007;56:183-197. <https://pubmed.ncbi.nlm.nih.gov/17881060/>

## 5. Neurological prognosis — the duration of coma decides it

> The basis for the `NEURO` cumulative injury integrator and the `IQEST` output. The model integrates injury
> as `(ammonia − 150) × time × (1 + 2×brain water)`.

32. Msall M, Batshaw ML, Suss R, Brusilow SW, Mellits ED. **Neurologic outcome in children with inborn errors of urea synthesis: outcome of urea-cycle enzymopathies.** *N Engl J Med* 1984;310:1500-1505. — the paper that first showed the **duration of coma** in neonatal hyperammonaemia to be the strongest predictor of later IQ. <https://pubmed.ncbi.nlm.nih.gov/6717538/>
33. Bachmann C. **Outcome and survival of 88 patients with urea cycle disorders: a retrospective evaluation.** *Eur J Pediatr* 2003;162:410-416. <https://pubmed.ncbi.nlm.nih.gov/12684900/>
34. Krivitzky L, et al. **Intellectual, adaptive, and behavioral functioning in children with urea cycle disorders.** *Pediatr Res* 2009;66:96-101. <https://pubmed.ncbi.nlm.nih.gov/19287347/>
35. Enns GM, et al. **Survival after treatment with phenylacetate and benzoate for urea-cycle disorders.** *N Engl J Med* 2007;356:2282-2292. — the Ammonul registry data; the relationship between duration of coma and survival. <https://pubmed.ncbi.nlm.nih.gov/17538087/>
36. Waisbren SE, et al. **Neuropsychological outcomes in individuals with urea cycle disorders.** *Mol Genet Metab* 2016;119:37-45. <https://pubmed.ncbi.nlm.nih.gov/27380995/>
37. Posset R, et al. **Impact of diagnosis and therapy on cognitive function in urea cycle disorders.** *Ann Neurol* 2019;86:116-128. <https://pubmed.ncbi.nlm.nih.gov/31018026/>

## 6. Nitrogen scavengers — the PK/PD of the phenylbutyrate axis

> Fixes `KHYDG` (GPB hydrolysis), `CLPBA`·`FMPAA` (β-oxidation), `VMPAGN`·`KMPAA`·`KMGLNP`
> (GLYATL1 conjugation), `CLPGN` (OAT secretion), and the difference between the NaPBA and GPB profiles.

38. Brusilow SW, Danney M, Waber LJ, et al. **Treatment of episodic hyperammonemia in children with inborn errors of urea synthesis.** *N Engl J Med* 1984;310:1630-1634. — the origin of benzoate and phenylacetate therapy. <https://pubmed.ncbi.nlm.nih.gov/6427608/>
39. Brusilow SW. **Phenylacetylglutamine may replace urea as a vehicle for waste nitrogen excretion.** *Pediatr Res* 1991;29:147-150. — the key stoichiometry that PAGN carries **two nitrogens** per molecule. <https://pubmed.ncbi.nlm.nih.gov/1903526/>
40. Diaz GA, Krivitzky LS, Mokhtarani M, et al. **Ammonia control and neurocognitive outcome among urea cycle disorder patients treated with glycerol phenylbutyrate.** *Hepatology* 2013;57:2171-2179. — the GPB vs NaPBA crossover design. **At the same molar dose of PBA, GPB gives a lower Cmax and a flatter 24-hour ammonia profile.** The anchor for scenarios 3 vs 4 of the model. <https://pubmed.ncbi.nlm.nih.gov/22961727/>
41. Smith W, Diaz GA, Lichter-Konecki U, et al. **Ammonia control in children ages 2 months through 5 years with urea cycle disorders: comparison of sodium phenylbutyrate and glycerol phenylbutyrate.** *J Pediatr* 2013;162:1228-1234. <https://pubmed.ncbi.nlm.nih.gov/23324524/>
42. Monteleone JPR, Mokhtarani M, Diaz GA, et al. **Population pharmacokinetic modeling and dosing simulations of nitrogen-scavenging compounds: disposition of glycerol phenylbutyrate and sodium phenylbutyrate in adult and pediatric patients with urea cycle disorders.** *J Clin Pharmacol* 2013;53:699-710. — the direct basis for the scavenger PK structure of this model (prodrug → PBA → PAA → PAGN) and for the U-PAGN recovery of ~60-75%. <https://pubmed.ncbi.nlm.nih.gov/23686462/>
43. Mokhtarani M, Diaz GA, Rhead W, et al. **Urinary phenylacetylglutamine as dosing biomarker for patients with urea cycle disorders.** *Mol Genet Metab* 2012;107:308-314. — the basis for using U-PAGN as a dose titration biomarker. <https://pubmed.ncbi.nlm.nih.gov/22958974/>
44. Mokhtarani M, Diaz GA, Rhead W, et al. **Elevated phenylacetic acid levels do not correlate with adverse events in patients with urea cycle disorders or hepatic encephalopathy and can be predicted based on the plasma PAA to PAGN ratio.** *Mol Genet Metab* 2013;110:446-453. — **a PAA:PAGN ratio > 2.5 (on a µg/mL basis) = saturation of conjugation**. The basis for the model's `PARATIO` output and for scenario 13. <https://pubmed.ncbi.nlm.nih.gov/24144944/>
45. Thibault A, et al. **Phase I study of phenylacetate administered twice daily to patients with cancer.** *Cancer* 1995;75:2932-2938. — the blood concentration range at which PAA neurotoxicity (somnolence, nausea) appears (≈500 µg/mL). <https://pubmed.ncbi.nlm.nih.gov/7773945/>
46. Berry SA, et al. **Safety and efficacy of glycerol phenylbutyrate for management of urea cycle disorders in patients aged 2 months to 2 years.** *Mol Genet Metab* 2017;122:46-53. <https://pubmed.ncbi.nlm.nih.gov/28916119/>
47. Scaglia F, Carter S, O'Brien WE, Lee B. **Effect of alternative pathway therapy on branched chain amino acid metabolism in urea cycle disorder patients.** *Mol Genet Metab* 2004;81(Suppl 1):S79-S85. — the phenylbutyrate-specific **BCAA depletion**. The `BCAAP` compartment and the `KBPAA` parameter of the model. <https://pubmed.ncbi.nlm.nih.gov/15050979/>
48. Burrage LC, et al. **Sodium phenylbutyrate decreases plasma branched-chain amino acids in patients with urea cycle disorders.** *Mol Genet Metab* 2014;113:131-135. <https://pubmed.ncbi.nlm.nih.gov/25042691/>

## 7. The benzoate-glycine axis, carglumic acid, cycle substrates

> Fixes `VMHIP`·`KMBZ`·`KMGLY` (GLYAT conjugation, **one nitrogen**), glycine depletion,
> the substitution of carglumic acid for NAGS, and the OTC-bypassing nitrogen excretion of citrulline (`KCITN`).

49. Batshaw ML, MacArthur RB, Tuchman M. **Alternative pathway therapy for urea cycle disorders: twenty years later.** *J Pediatr* 2001;138(1 Suppl):S46-S55. — a comparison of the nitrogen stoichiometry of benzoate against that of phenylbutyrate. <https://pubmed.ncbi.nlm.nih.gov/11148549/>
50. Kasapkara ÇS, et al. **N-carbamylglutamate treatment for acute neonatal hyperammonemia in isolated methylmalonic acidemia.** *Eur J Pediatr* 2011;170:799-801. <https://pubmed.ncbi.nlm.nih.gov/21170724/>
51. Daniotti M, la Marca G, Fiorini P, Filippi L. **New developments in the treatment of hyperammonemia: emerging use of carglumic acid.** *Int J Gen Med* 2011;4:21-28. <https://pubmed.ncbi.nlm.nih.gov/21403788/>
52. Caldovic L, et al. **N-acetylglutamate synthase deficiency: a heterogeneous disorder.** *Mol Genet Metab* 2007;91:36-41. — establishes NAGS deficiency as the only UCD that carglumic acid **phenotypically cures**. Scenario 11 of the model. <https://pubmed.ncbi.nlm.nih.gov/17368065/>
53. Brusilow SW. **Arginine, an indispensable amino acid for patients with inborn errors of urea synthesis.** *J Clin Invest* 1984;74:2144-2148. <https://pubmed.ncbi.nlm.nih.gov/6511918/>
54. Nagasaka H, et al. **Effects of arginine treatment on nutrition, growth and urea cycle function in seven Japanese boys with late-onset ornithine transcarbamylase deficiency.** *Eur J Pediatr* 2006;165:618-624. <https://pubmed.ncbi.nlm.nih.gov/16691408/>
55. Erez A, Nagamani SCS, Shchelochkov OA, et al. **Requirement of argininosuccinate lyase for systemic nitric oxide production.** *Nat Med* 2011;17:1619-1626. — the discovery that ASL is an obligatory component of the NOS complex. The `NO_deficit` and `Hypertension_ASL` nodes of the model. <https://pubmed.ncbi.nlm.nih.gov/22081021/>
56. Nagamani SCS, et al. **Nitric-oxide supplementation for treatment of long-term complications in argininosuccinic aciduria.** *Am J Hum Genet* 2012;90:836-846. <https://pubmed.ncbi.nlm.nih.gov/22541557/>

## 8. Glutamine — both the buffer store and the early warning

> The `ucd_vgs()` function (a high-affinity plus a low-affinity GS component), the Gln:NH3 relationship,
> and the basis for the clinical observation that glutamine rises **before** ammonia.

57. Darmaun D, Matthews DE, Bier DM. **Glutamine and glutamate kinetics in humans.** *Am J Physiol* 1986;251:E117-E126. — whole-body glutamine flux ≈ 300 µmol/kg/h. The direct basis for the `KGSKG` and `VG0KG` calibration. <https://pubmed.ncbi.nlm.nih.gov/2873746/>
58. Lee B, Diaz GA, Rhead W, et al. **Blood ammonia and glutamine as predictors of hyperammonemic crises in patients with urea cycle disorder.** *Genet Med* 2015;17:561-568. — **a glutamine > 1000 µmol/L predicts a crisis**. The early warning logic of the model. <https://pubmed.ncbi.nlm.nih.gov/25341094/>
59. Maestri NE, McGowan KD, Brusilow SW. **Plasma glutamine concentration: a guide in the management of urea cycle disorders.** *J Pediatr* 1992;121:259-261. <https://pubmed.ncbi.nlm.nih.gov/1640295/>
60. Häussinger D. **Nitrogen metabolism in liver: structural and functional organization and physiological relevance.** *Biochem J* 1990;267:281-290. — periportal (high-capacity ureagenesis) / perivenous (high-affinity glutamine synthesis) compartmentalisation. The `Hepatocyte_periportal` / `Hepatocyte_perivenous` nodes of the model. <https://pubmed.ncbi.nlm.nih.gov/1970241/>
61. Weiner ID, Verlander JW. **Renal ammonia metabolism and transport.** *Compr Physiol* 2013;3:201-220. — excretion of ammonia derived from renal glutaminase (the `CLRGLN` pathway of the model). <https://pubmed.ncbi.nlm.nih.gov/23720285/>

## 9. Orotic acid — the signal that separates OTC from CPS1/NAGS

62. Bachmann C, Colombo JP. **Diagnostic value of orotic acid excretion in heritable disorders of the urea cycle and in hyperammonemia due to organic acidurias.** *Eur J Pediatr* 1980;134:109-113. <https://pubmed.ncbi.nlm.nih.gov/7449803/>
63. Salerno C, Crifò C. **Diagnostic value of urinary orotic acid levels: applicable separation methods.** *J Chromatogr B* 2002;781:57-71. <https://pubmed.ncbi.nlm.nih.gov/12450653/>

## 10. Emergency management · extracorporeal removal · rebound

> The basis for `CLHD` (dialysis clearance), `FHDGLN`, and the **post-dialysis rebound** that comes out of the
> deep glutamine store (`GLNM`). Scenario 9 vs 10.

64. Picca S, et al. **Extracorporeal dialysis in neonatal hyperammonemia: modalities and prognostic indicators.** *Pediatr Nephrol* 2001;16:862-867. <https://pubmed.ncbi.nlm.nih.gov/11685590/>
65. Spinale JM, et al. **High-dose continuous renal replacement therapy for neonatal hyperammonemia.** *Pediatr Nephrol* 2013;28:983-986. <https://pubmed.ncbi.nlm.nih.gov/23515666/>
66. Schaefer F, et al. **Dialysis in neonates with inborn errors of metabolism.** *Nephrol Dial Transplant* 1999;14:910-918. <https://pubmed.ncbi.nlm.nih.gov/10328466/>
67. Wiegand C, Thompson T, Bock GH, Mathis RK, Kjellstrand CM, Mauer SM. **The management of life-threatening hyperammonemia: a comparison of several therapeutic modalities.** *J Pediatr* 1980;96:142-144. — haemodialysis is overwhelmingly faster than peritoneal dialysis or exchange transfusion. <https://pubmed.ncbi.nlm.nih.gov/7350297/>
68. Ah Mew N, et al. **Comparison of sodium phenylbutyrate and glycerol phenylbutyrate in the acute management of hyperammonemia.** *Mol Genet Metab* 2018;124:1-6. <https://pubmed.ncbi.nlm.nih.gov/29655841/>

## 11. Liver transplantation · gene therapy · mRNA therapy

> The basis for `ACTX`/`TTX` (transplantation) and `TRANSG`·`OTCX`·`KEXPR` (AAV/ mRNA).
> Scenario 12.

69. Yu L, et al. **Liver transplantation for urea cycle disorders: analysis of the United Network for Organ Sharing database.** *Transplant Proc* 2015;47:2413-2418. <https://pubmed.ncbi.nlm.nih.gov/26518941/>
70. Kido J, et al. **Long-term outcome and intervention of urea cycle disorders in Japan.** *J Inherit Metab Dis* 2012;35:777-785. — transplantation phenotypically cures ureagenesis but **does not reverse existing CNS injury**. <https://pubmed.ncbi.nlm.nih.gov/22167275/>
71. Raper SE, et al. **Fatal systemic inflammatory response syndrome in a ornithine transcarbamylase deficient patient following adenoviral gene transfer.** *Mol Genet Metab* 2003;80:148-158. — the turning point in the history of gene therapy (the Gelsinger case). The `Immune_AAV` node of the model. <https://pubmed.ncbi.nlm.nih.gov/14567964/>
72. Wang L, et al. **AAV gene therapy corrects OTC deficiency and prevents liver fibrosis in aged neonatal ornithine transcarbamylase-deficient mice.** *Mol Genet Metab* 2017;120:299-305. <https://pubmed.ncbi.nlm.nih.gov/28202336/>
73. Harding CO, et al. **Safety and efficacy of DTX301, an AAV8-mediated gene transfer for adults with late-onset ornithine transcarbamylase deficiency: interim results.** *Mol Genet Metab* 2021;132:S49. — recovery of ureagenesis flux and reduction of the scavenger dose. (Programme overview) <https://pubmed.ncbi.nlm.nih.gov/33642232/>
74. Prieve MG, et al. **Targeted mRNA therapy for ornithine transcarbamylase deficiency.** *Mol Ther* 2018;26:801-813. — repeated LNP-mRNA dosing restores hepatic OTC activity in pulses. The mRNA mode of the model. <https://pubmed.ncbi.nlm.nih.gov/29433939/>
75. Diez-Fernandez C, Häberle J. **Targeting CPS1 in the treatment of carbamoyl phosphate synthetase 1 deficiency.** *Expert Opin Ther Targets* 2017;21:391-399. <https://pubmed.ncbi.nlm.nih.gov/28281899/>

## 12. Precipitants · drug interactions (valproate and others)

> Scenario 14 and the `Valproate`, `Corticosteroid_DDI`, and `Probenecid` nodes of the map.

76. Coulter DL, Allen RJ. **Secondary hyperammonaemia: a possible mechanism for valproate encephalopathy.** *Lancet* 1980;1:1310-1311. <https://pubmed.ncbi.nlm.nih.gov/6104119/>
77. Aires CCP, et al. **New insights on the mechanisms of valproate-induced hyperammonemia: inhibition of hepatic N-acetylglutamate synthase activity by valproyl-CoA.** *J Hepatol* 2011;55:426-434. — valproate directly inhibits NAGS. The `VPAI` parameter of the model. <https://pubmed.ncbi.nlm.nih.gov/21147182/>
78. Tuchman M, Yudkoff M. **Blood levels of ammonia and nitrogen scavenging amino acids in patients with inherited hyperammonemia.** *Mol Genet Metab* 1999;66:10-15. <https://pubmed.ncbi.nlm.nih.gov/9973543/>
79. Nott L, Price TJ, Pittman K, Patterson K, Fletcher J. **Hyperammonemia encephalopathy: an important cause of neurological deterioration following chemotherapy.** *Leuk Lymphoma* 2007;48:1702-1711. <https://pubmed.ncbi.nlm.nih.gov/17786705/>
80. Lipskind S, Loanzon S, Simi E, Ouyang DW. **Hyperammonemic coma in ornithine transcarbamylase deficiency: a case of postpartum decompensation.** *Obstet Gynecol* 2011;117:503-505. — postpartum uterine involution as a representative catabolic precipitant. <https://pubmed.ncbi.nlm.nih.gov/21252804/>

## 13. Nutrition · protein prescription · growth

> The basis for the `PROT`, `FDISP`, and `NOBLKG` parameters and for the `PROTTOL` (natural protein tolerance) output.

81. Adam S, et al. **Dietary management of urea cycle disorders: European practice.** *Mol Genet Metab* 2013;110:439-445. <https://pubmed.ncbi.nlm.nih.gov/24113687/>
82. Singh RH. **Nutritional management of patients with urea cycle disorders.** *J Inherit Metab Dis* 2007;30:880-887. <https://pubmed.ncbi.nlm.nih.gov/17957501/>
83. Boyle M, et al. **Growth in patients with urea cycle disorders.** *Mol Genet Metab* 2014;113:220-224. <https://pubmed.ncbi.nlm.nih.gov/25266922/>
84. WHO/FAO/UNU. **Protein and amino acid requirements in human nutrition.** *WHO Tech Rep Ser* 935, 2007. — the source of the obligatory nitrogen loss value. <https://pubmed.ncbi.nlm.nih.gov/18330140/>

## 14. QSP · mrgsolve methodology

85. Baron KT, et al. **mrgsolve: Simulate from ODE-Based Models.** R package. <https://mrgsolve.org/>
86. Nijsen MJMA, et al. **Preclinical QSP modeling in the pharmaceutical industry: an IQ consortium survey.** *CPT Pharmacometrics Syst Pharmacol* 2018;7:135-146. <https://pubmed.ncbi.nlm.nih.gov/29349875/>
87. Bai JPF, et al. **Quantitative systems pharmacology: landscape analysis of regulatory submissions to the US FDA.** *CPT Pharmacometrics Syst Pharmacol* 2021;10:1479-1484. <https://pubmed.ncbi.nlm.nih.gov/34617412/>
88. Gadkar K, et al. **A six-stage workflow for robust application of systems pharmacology.** *CPT Pharmacometrics Syst Pharmacol* 2016;5:235-249. <https://pubmed.ncbi.nlm.nih.gov/27299936/>

---

## Appendix A — Which reference each model parameter came from

| Parameter | Value | Basis |
|----------|-----|------|
| `UCAP` (maximum ureagenesis capacity) | 2400 µmol N/kg/h | the maximum ureagenesis rate is 4-5 times the habitual nitrogen load [9, 13, 15] |
| `KMNH4` | 120 µmol/L | an internally consistent value, making a healthy person use ~20% of capacity at an ammonia of 30 µmol/L [9, 29] |
| `VG0KG` / `KGSKG` / `KGSAT` | 184.6 / 3.843 / 1200 | whole-body glutamine flux ≈300 µmol/kg/h [57], the observed Gln-NH3 relationship [58, 59] |
| `CLRGKG` | 0.0275 L/h/kg | excretion of renal glutamine-derived nitrogen of 30-50 mmol/day [61] |
| `VMPKG` / `KMPAA` / `KMGLNP` | 115 µmol/kg/h / 300 / 400 | U-PAGN recovery of 60-75% [42, 43] |
| `KHYDG` (GPB hydrolysis) | 0.30 /h | the low Cmax and flat profile of GPB [40, 42] |
| PAA toxicity threshold | 500 µg/mL | [44, 45] |
| PAA:PAGN threshold | 2.5 (µg/mL) | [44] |
| `PSBBB` / `PKANH3` / `PHBR` | 400 L/h / 9.15 / 7.05 | non-ionic diffusion plus ion trapping in the brain [29] |
| `VMGSB` / `KOUTGB` | 892 µmol/h / 0.05 /h | a rise in brain glutamine from 5 to 15-20 mmol/L [25, 26, 27] |
| `KINMI` / `KEXTMI` / `FCOMP` | 0.10 / 0.02 / 0.60 | the **partial** osmotic compensation by myo-inositol [25, 22] |
| `KINJ` / `IQ50` | 0.10 / 50 | the duration-of-coma to IQ relationship [32, 33, 37] |
| Ammonia thresholds 100 / 200 / 360 | — | [1] |
| `NOBLKG` (obligatory nitrogen outflow) | 43 µmol N/kg/h | obligatory nitrogen loss on a protein-free diet [84] |

## Appendix B — Clinical observations the model reproduces

| Observation | How it is reproduced |
|------|-----------|
| The disease is not a slope but a **cliff** | Michaelis-Menten ureagenesis + the glutamine buffer store → ammonia 484 at 5% residual activity, 61 µmol/L at 40% (the cliff plot on tab 1) |
| Glutamine rises **before** ammonia | the buffering structure by which `ucd_vgs()` moves nitrogen out of ammonia into glutamine [58] |
| A patient with **chronic** hyperammonaemia is more vulnerable to an acute precipitant | `MINSB` already depleted at baseline → no osmotic buffer [25] |
| GPB has a **lower Cmax** at the same molar dose | the lipase-dependent `GPBG → PBAG` hydrolysis step [40] |
| Benzoate carries one nitrogen per mole, phenylbutyrate two | `-vhip` (1 N) in `dxdt_NH4C` vs the PAGN pathway (2 N) [39, 49] |
| In a neonatal crisis **dialysis is more decisive than drugs** | scenario 9 vs 10 (ammonia 1182 → 236 vs 1182 → 31) [67] |
| **Post-dialysis rebound** | redistribution from the deep glutamine store `GLNM` (scenario 10, 86-120 h) |
| Carglumic acid **phenotypically cures** NAGS deficiency | `GNCG` normalises `f_NAG` → protein tolerance 0.46 → 1.35 g/kg/d (scenario 11) [52] |
| What decides prognosis is not the mode of treatment but the **delay to diagnosis** | in scenarios 9/10, after a 60-hour delay even dialysis recovers only 3 IQ points [32] |

---

> **Disclaimer.** This model and this collation of the literature are for educational and research purposes.
> The parameters are approximations derived from the published literature and have not undergone independent verification.
> They must not be used for clinical decision-making, for prescribing, or for regulatory submission.
