# Chronic hyperkalaemia QSP model — references

**Chronic hyperkalaemia in CKD and heart failure — the RAAS-inhibitor potassium dilemma**

This document collects, for `hk_qsp_model.dot`, `hk_mrgsolve_model.R`, `hk_reference_model.py`, and
`hk_shiny_app_en.R`, **the provenance of every structural assumption and parameter** used in them.
Each entry also indicates which part of the model it was used in.

> Notation: **[F]** = a parameter was *fitted* to this reference · **[V]** = used for *validation*
> only, not fitted · **[S]** = the basis of a *structure* ·
> **[C]** = clinical *context*

---

## 1. The basic structure of potassium homeostasis — intra- and extracellular distribution and buffering

The basis for the core relation of the model,
`K_total = Ce·V_ECF + Ci0·LAMrel·(Ce/Ce0)^α·V_ICF`, and for α = 0.25.

1. **[S]** Palmer BF. *Regulation of Potassium Homeostasis.* Clin J Am Soc Nephrol. 2015;10(6):1050-60. — https://pubmed.ncbi.nlm.nih.gov/24721891/
2. **[S]** Palmer BF, Clegg DJ. *Physiology and Pathophysiology of Potassium Homeostasis: Core Curriculum 2019.* Am J Kidney Dis. 2019;74(5):682-95. — https://pubmed.ncbi.nlm.nih.gov/31227226/
3. **[S]** Youn JH, McDonough AA. *Recent advances in understanding integrative control of potassium homeostasis.* Annu Rev Physiol. 2009;71:381-401. — https://pubmed.ncbi.nlm.nih.gov/18759636/
4. **[V]** Sterns RH, Cox M, Feig PU, Singer I. *Internal potassium balance and the control of the plasma potassium concentration.* Medicine (Baltimore). 1981;60(5):339-54. — https://pubmed.ncbi.nlm.nih.gov/6268928/ — a total body deficit of 200-400 mmol at a serum K of 3.0. The model **does not look at** this value and predicts 302 mmol (derived from α).
5. **[S]** Clausen T. *Quantification of Na+,K+ pumps and their transport rate in skeletal muscle: functional significance.* J Gen Physiol. 2013;142(4):327-45. — https://pubmed.ncbi.nlm.nih.gov/24081980/ — skeletal muscle holds ~75% of total body K and is the seat of the buffering.
6. **[S]** McDonough AA, Youn JH. *Potassium Homeostasis: The Knowns, the Unknowns, and the Health Benefits.* Physiology. 2017;32(2):100-11. — https://pubmed.ncbi.nlm.nih.gov/28202621/
7. **[S]** Gumz ML, Rabinowitz L, Wingo CS. *An Integrated View of Potassium Homeostasis.* N Engl J Med. 2015;373(1):60-72. — https://pubmed.ncbi.nlm.nih.gov/26132942/

## 2. Regulation of transcellular shift — insulin · β2 · pH · osmolality

Each term of `LAMrel = f_ins · f_β2 · f_pH · f_glu · f_aldo`.

8. **[F]** Nguyen TQ, Maalouf NM, Sakhaee K, Moe OW. *Comparison of insulin action on glucose versus potassium uptake in humans.* Clin J Am Soc Nephrol. 2011;6(7):1533-9. — https://pubmed.ncbi.nlm.nih.gov/21734082/ — insulin-mediated K uptake saturates at a lower concentration than glucose uptake → basal insulin already does most of the work (the basis for making `f_ins` a function over the whole range in the model).
9. **[V]** Allon M, Copkney C. *Albuterol and insulin for treatment of hyperkalemia in hemodialysis patients.* Kidney Int. 1990;38(5):869-72. — https://pubmed.ncbi.nlm.nih.gov/2266671/ — insulin -0.6 to -1.0, salbutamol -0.6 to -1.0 mmol/L.
10. **[V]** Ngugi NN, McLigeyo SO, Kayima JK. *Treatment of hyperkalaemia by altering the transcellular gradient in patients with renal failure.* East Afr Med J. 1997;74(8):503-9. — https://pubmed.ncbi.nlm.nih.gov/9487416/
11. **[S]** Aronson PS, Giebisch G. *Effects of pH on potassium: new explanations for old observations.* J Am Soc Nephrol. 2011;22(11):1981-9. — https://pubmed.ncbi.nlm.nih.gov/21980112/ — that only mineral acidosis shifts K appreciably, and organic acidosis (lactic acid, ketoacids) does not.
12. **[S]** Adrogué HJ, Madias NE. *Changes in plasma potassium concentration during acute acid-base disturbances.* Am J Med. 1981;71(3):456-67. — https://pubmed.ncbi.nlm.nih.gov/7025622/
13. **[S]** Clausen T. *Hormonal and pharmacological modification of plasma potassium homeostasis.* Fundam Clin Pharmacol. 2010;24(5):595-605. — https://pubmed.ncbi.nlm.nih.gov/20618871/
14. **[C]** Kamel KS, Wei C. *Controversial issues in the treatment of hyperkalaemia.* Nephrol Dial Transplant. 2003;18(11):2215-8. — https://pubmed.ncbi.nlm.nih.gov/14551344/

## 3. Renal potassium handling and the distal tubule (the ASDN)

The structure of `E_ren = FD·filt + S_cap·RASDN·(Ce/(Km+Ce))·q^0.5·(HCO3/24)^n`.

15. **[S]** Palmer BF, Clegg DJ. *Hyperkalemia across the Continuum of Kidney Function.* Clin J Am Soc Nephrol. 2018;13(1):155-7. — https://pubmed.ncbi.nlm.nih.gov/29114006/
16. **[F]** Hayes CP Jr, McLeod ME, Robinson RR. *An extravenal [sic] mechanism for the maintenance of potassium balance in severe chronic renal failure.* Trans Assoc Am Physicians. 1967;80:207-16. — https://pubmed.ncbi.nlm.nih.gov/6082243/ — upregulation of colonic secretion in renal failure (the model's `KC_COL`).
17. **[F]** Schultze RG, Taggart DD, Shapiro H, et al. *On the adaptation in potassium excretion associated with nephron reduction in the dog.* J Clin Invest. 1971;50(5):1061-8. — https://pubmed.ncbi.nlm.nih.gov/5552407/ — upregulation of secretory capacity per remaining nephron (the model's `ADAPT_P` and `ADAPT_MX`).
18. **[S]** Welling PA. *Roles and Regulation of Renal K Channels.* Annu Rev Physiol. 2016;78:415-35. — https://pubmed.ncbi.nlm.nih.gov/26654186/ — the division of labour between ROMK and BK (Maxi-K).
19. **[S]** Wang WH, Giebisch G. *Regulation of potassium (K) handling in the renal collecting duct.* Pflugers Arch. 2009;458(1):157-68. — https://pubmed.ncbi.nlm.nih.gov/18839206/
20. **[S]** Terker AS, Zhang C, McCormick JA, et al. *Potassium modulates electrolyte balance and blood pressure through effects on distal cell voltage and chloride.* Cell Metab. 2015;21(1):39-50. — https://pubmed.ncbi.nlm.nih.gov/25565204/ — the Kir4.1/5.1–WNK–SPAK–NCC switch (cluster 5 of the map).
21. **[S]** Cuevas CA, Su XT, Wang MX, et al. *Potassium Sensing by Renal Distal Tubules Requires Kir4.1.* J Am Soc Nephrol. 2017;28(6):1814-25. — https://pubmed.ncbi.nlm.nih.gov/28052988/
22. **[S]** Sorensen MV, Grossmann S, Roesinger M, et al. *Rapid dephosphorylation of the renal sodium chloride cotransporter in response to oral potassium intake in mice.* Kidney Int. 2013;83(5):811-24. — https://pubmed.ncbi.nlm.nih.gov/23447069/
23. **[S]** Preston RA, Afshartous D, Rodco R, et al. *Evidence for a gastrointestinal-renal kaliuretic signaling axis in humans.* Kidney Int. 2015;88(6):1383-91. — https://pubmed.ncbi.nlm.nih.gov/26308672/ — the gut-kidney feedback (the map's `gutsensor`).
24. **[S]** Rabelink TJ, Koomans HA, Hené RJ, Dorhout Mees EJ. *Early and late adjustment to potassium loading in humans.* Kidney Int. 1990;38(5):942-7. — https://pubmed.ncbi.nlm.nih.gov/2266680/
25. **[S]** DuBose TD Jr. *Regulation of Potassium Homeostasis in CKD.* Adv Chronic Kidney Dis. 2017;24(5):305-14. — https://pubmed.ncbi.nlm.nih.gov/29031357/

## 4. Acid-base and potassium

The basis for the model's conclusion that **alkali therapy acts through renal secretion rather than through cellular shift**.

26. **[F]** de Brito-Ashurst I, Varagunam M, Raftery MJ, Yaqoob MM. *Bicarbonate supplementation slows progression of CKD and improves nutritional status.* J Am Soc Nephrol. 2009;20(9):2075-84. — https://pubmed.ncbi.nlm.nih.gov/19608703/
27. **[F]** Di Iorio BR, Bellasi A, Raphael KL, et al. (UBI Study) *Treatment of metabolic acidosis with sodium bicarbonate delays progression of chronic kidney disease.* J Nephrol. 2019;32(6):989-1001. — https://pubmed.ncbi.nlm.nih.gov/31598912/
28. **[S]** Wesson DE, Buysse JM, Bushinsky DA. *Mechanisms of Metabolic Acidosis-Induced Kidney Injury in CKD.* J Am Soc Nephrol. 2020;31(3):469-82. — https://pubmed.ncbi.nlm.nih.gov/31988269/
29. **[V]** Blumberg A, Weidmann P, Shaw S, Gnädinger M. *Effect of various therapeutic approaches on plasma potassium and major regulating factors in terminal renal failure.* Am J Med. 1988;85(4):507-12. — https://pubmed.ncbi.nlm.nih.gov/3052050/ — the classic negative result that **bicarbonate is not a treatment for acute hyperkalaemia**.
30. **[V]** Allon M, Shanklin N. *Effect of bicarbonate administration on plasma potassium in dialysis patients.* Am J Kidney Dis. 1996;28(4):508-14. — https://pubmed.ncbi.nlm.nih.gov/8840939/

## 5. The RAAS axis and the mineralocorticoid receptor

31. **[S]** Arroyo JP, Ronzaud C, Lagnaz D, et al. *Aldosterone paradox: differential regulation of ion transport in distal nephron.* Physiology. 2011;26(2):115-23. — https://pubmed.ncbi.nlm.nih.gov/21487030/ — the structure of the aldosterone paradox (Na conservation vs K secretion).
32. **[S]** Himathongkam T, Dluhy RG, Williams GH. *Potassim[sic]-aldosterone-renin interrelationships.* J Clin Endocrinol Metab. 1975;41(1):153-9. — https://pubmed.ncbi.nlm.nih.gov/1167307/ — direct stimulation of aldosterone by serum K (the model's `SK_ALDO`, the aldosterone escape pathway).
33. **[C]** DeFronzo RA. *Hyperkalemia and hyporeninemic hypoaldosteronism.* Kidney Int. 1980;17(1):118-34. — https://pubmed.ncbi.nlm.nih.gov/6990088/ — type 4 renal tubular acidosis.
34. **[F]** Pitt B, Zannad F, Remme WJ, et al. (RALES) *The effect of spironolactone on morbidity and mortality in patients with severe heart failure.* N Engl J Med. 1999;341(10):709-17. — https://pubmed.ncbi.nlm.nih.gov/10471456/ — a K rise of +0.30 mmol/L on spironolactone 25 mg (the fitting anchor for the model's `KI_MRA`).
35. **[V]** Pitt B, Remme W, Zannad F, et al. (EPHESUS) *Eplerenone, a selective aldosterone blocker, in patients with left ventricular dysfunction after myocardial infarction.* N Engl J Med. 2003;348(14):1309-21. — https://pubmed.ncbi.nlm.nih.gov/12668699/
36. **[V]** Bakris GL, Agarwal R, Anker SD, et al. (FIDELIO-DKD) *Effect of Finerenone on Chronic Kidney Disease Outcomes in Type 2 Diabetes.* N Engl J Med. 2020;383(23):2219-29. — https://pubmed.ncbi.nlm.nih.gov/33264825/ — an observed K rise of +0.23; the model gives +0.10 at the assumed MR load, a **mismatch**, and reports it as such.
37. **[V]** Pitt B, Filippatos G, Agarwal R, et al. (FIGARO-DKD) *Cardiovascular Events with Finerenone in Kidney Disease and Type 2 Diabetes.* N Engl J Med. 2021;385(24):2252-63. — https://pubmed.ncbi.nlm.nih.gov/34449181/
38. **[C]** Juurlink DN, Mamdani MM, Lee DS, et al. *Rates of hyperkalemia after publication of the Randomized Aldactone Evaluation Study.* N Engl J Med. 2004;351(6):543-51. — https://pubmed.ncbi.nlm.nih.gov/15295047/ — the surge in hyperkalaemia in the population outside the trial: a real instance of the "the eGFR band the trial excluded" problem the model emphasises.

## 6. Potassium binders

39. **[V]** Weir MR, Bakris GL, Bushinsky DA, et al. (OPAL-HK) *Patiromer in patients with kidney disease and hyperkalemia receiving RAAS inhibitors.* N Engl J Med. 2015;372(3):211-21. — https://pubmed.ncbi.nlm.nih.gov/25415805/ — ΔK -1.01 at 4 weeks (model -0.95, not fitted).
40. **[V]** Bakris GL, Pitt B, Weir MR, et al. (AMETHYST-DN) *Effect of Patiromer on Serum Potassium Level in Patients With Hyperkalemia and Diabetic Kidney Disease.* JAMA. 2015;314(2):151-61. — https://pubmed.ncbi.nlm.nih.gov/26172895/
41. **[V]** Kosiborod M, Rasmussen HS, Lavin P, et al. (HARMONIZE) *Effect of sodium zirconium cyclosilicate on potassium lowering for 28 days among outpatients with hyperkalemia.* JAMA. 2014;312(21):2223-33. — https://pubmed.ncbi.nlm.nih.gov/25402495/ — K of 4.8/4.5/4.4 on maintenance doses of 5/10/15 g (model 4.89/4.63/4.49).
42. **[V]** Packham DK, Rasmussen HS, Lavin PT, et al. (ZS-003) *Sodium zirconium cyclosilicate in hyperkalemia.* N Engl J Med. 2015;372(3):222-31. — https://pubmed.ncbi.nlm.nih.gov/25415807/
43. **[S]** Stavros F, Yang A, Leon A, et al. *Characterization of structure and function of ZS-9, a K selective ion trap.* PLoS One. 2014;9(12):e114686. — https://pubmed.ncbi.nlm.nih.gov/25531770/ — the structural basis for SZC acting in the stomach and proximal gut (why its site of action is set differently from patiromer's in the model).
44. **[S]** Li L, Harrison SD, Cope MJ, Park C, et al. *Mechanism of Action and Pharmacology of Patiromer, a Nonabsorbed Cross-Linked Polymer That Lowers Serum Potassium Concentration.* J Cardiovasc Pharmacol Ther. 2016;21(5):456-65. — https://pubmed.ncbi.nlm.nih.gov/26856345/ — Ca-K exchange in the distal colon.
45. **[C]** Agarwal R, Rossignol P, Romero A, et al. (AMBER) *Patiromer versus placebo to enable spironolactone use in patients with resistant hypertension and CKD.* Lancet. 2019;394(10208):1540-50. — https://pubmed.ncbi.nlm.nih.gov/31533906/ — the randomised basis for the concept of "RAASi enablement".
46. **[C]** Butler J, Anker SD, Lund LH, et al. (DIAMOND) *Patiromer for the management of hyperkalaemia in heart failure with reduced ejection fraction.* Eur Heart J. 2022;43(41):4362-73. — https://pubmed.ncbi.nlm.nih.gov/35900838/
47. **[C]** Sterns RH, Rojas M, Bernstein P, Chennupati S. *Ion-exchange resins for the treatment of hyperkalemia: are they safe and effective?* J Am Soc Nephrol. 2010;21(5):733-5. — https://pubmed.ncbi.nlm.nih.gov/20167700/ — the weakness of the evidence for SPS (polystyrene sulfonate) and the signal of colonic necrosis.
48. **[C]** Natale P, Palmer SC, Ruospo M, et al. *Potassium binders for chronic hyperkalaemia in people with CKD.* Cochrane Database Syst Rev. 2020;6:CD013165. — https://pubmed.ncbi.nlm.nih.gov/32588430/

## 7. Acute management of hyperkalaemia

49. **[C]** Lindner G, Burdmann EA, Clase CM, et al. *Acute hyperkalemia in the emergency department: a summary from a Kidney Disease: Improving Global Outcomes conference.* Eur J Emerg Med. 2020;27(5):329-37. — https://pubmed.ncbi.nlm.nih.gov/32852924/
50. **[V]** Harel Z, Kamel KS. *Optimal Dose and Method of Administration of Intravenous Insulin in the Management of Emergency Hyperkalemia: A Systematic Review.* PLoS One. 2016;11(5):e0154963. — https://pubmed.ncbi.nlm.nih.gov/27148740/
51. **[C]** Coca A, Valencia AL, Bustamante J, et al. *Hypoglycemia following intravenous insulin plus glucose for hyperkalemia in patients with impaired renal function.* PLoS One. 2017;12(2):e0172961. — https://pubmed.ncbi.nlm.nih.gov/28245289/ — the 15-20% incidence of hypoglycaemia, which the model **does not reproduce** (stated as a limitation).
52. **[S]** Parham WA, Mehdirad AA, Biermann KM, Fredman CS. *Hyperkalemia revisited.* Tex Heart Inst J. 2006;33(1):40-7. — https://pubmed.ncbi.nlm.nih.gov/16572868/
53. **[S]** Diercks DB, Shumaik GM, Harrigan RA, et al. *Electrocardiographic manifestations: electrolyte abnormalities.* J Emerg Med. 2004;27(2):153-60. — https://pubmed.ncbi.nlm.nih.gov/15261358/
54. **[V]** Montague BT, Ouellette JR, Buller GK. *Retrospective review of the frequency of ECG changes in hyperkalemia.* Clin J Am Soc Nephrol. 2008;3(2):324-30. — https://pubmed.ncbi.nlm.nih.gov/18235147/ — the low sensitivity of peaked T waves (why the model does not use the T wave as a diagnostic index).
55. **[S]** Weiss JN, Qu Z, Shivkumar K. *Electrophysiology of Hypokalemia and Hyperkalemia.* Circ Arrhythm Electrophysiol. 2017;10(3):e004667. — https://pubmed.ncbi.nlm.nih.gov/28314851/ — the quantitative relationship between Em, Na channel availability, and conduction velocity (cluster 12 of the model).

## 8. Epidemiology and outcomes

56. **[F]** Kovesdy CP, Matsushita K, Sang Y, et al. (CKD Prognosis Consortium) *Serum potassium and adverse outcomes across the range of kidney function: a CKD Prognosis Consortium meta-analysis.* Eur Heart J. 2018;39(17):1535-42. — https://pubmed.ncbi.nlm.nih.gov/29554312/ — the U-shaped risk curve (the basis for fitting the model's `BK_HI` and `BK_LO`).
57. **[F]** Collins AJ, Pitt B, Reaven N, et al. *Association of Serum Potassium with All-Cause Mortality in Patients with and without Heart Failure, CKD, and/or Diabetes.* Am J Nephrol. 2017;46(3):213-21. — https://pubmed.ncbi.nlm.nih.gov/28866674/
58. **[C]** Einhorn LM, Zhan M, Hsu VD, et al. *The frequency of hyperkalemia and its significance in chronic kidney disease.* Arch Intern Med. 2009;169(12):1156-62. — https://pubmed.ncbi.nlm.nih.gov/19546417/
59. **[C]** Luo J, Brunelli SM, Jensen DE, Yang A. *Association between Serum Potassium and Outcomes in Patients with Reduced Kidney Function.* Clin J Am Soc Nephrol. 2016;11(1):90-100. — https://pubmed.ncbi.nlm.nih.gov/26500246/
60. **[C]** James G, Kim J, Mellström C, Ford KL, et al. *Serum potassium variability as a predictor of clinical outcomes in patients with cardiorenal disease or diabetes.* Clin Kidney J. 2022;15(4):758-70. — https://pubmed.ncbi.nlm.nih.gov/35371436/

## 9. The cost of RAASi down-titration and discontinuation — the central result of the model

61. **[F]** Epstein M, Reaven NL, Funk SE, et al. *Evaluation of the treatment gap between clinical guidelines and the utilization of renin-angiotensin-aldosterone system inhibitors.* Am J Manag Care. 2015;21(11 Suppl):S212-20. — https://pubmed.ncbi.nlm.nih.gov/26619183/ — the association between RAASi down-titration or discontinuation after hyperkalaemia and increased mortality.
62. **[V]** Xie X, Liu Y, Perkovic V, et al. *Renin-Angiotensin System Inhibitors and Kidney and Cardiovascular Outcomes in Patients With CKD: A Bayesian Network Meta-analysis.* Am J Kidney Dis. 2016;67(5):728-41. — https://pubmed.ncbi.nlm.nih.gov/26597926/ — the model's `LNHR_RD` and the eGFR slope effect.
63. **[V]** Brenner BM, Cooper ME, de Zeeuw D, et al. (RENAAL) *Effects of losartan on renal and cardiovascular outcomes in patients with type 2 diabetes and nephropathy.* N Engl J Med. 2001;345(12):861-9. — https://pubmed.ncbi.nlm.nih.gov/11565518/
64. **[V]** Fu EL, Evans M, Clase CM, et al. *Stopping Renin-Angiotensin System Inhibitors in Patients with Advanced CKD and Risk of Adverse Outcomes.* J Am Soc Nephrol. 2021;32(2):424-35. — https://pubmed.ncbi.nlm.nih.gov/33372009/
65. **[C]** Bhandari S, Mehta S, Khwaja A, et al. (STOP-ACEi) *Renin-Angiotensin System Inhibition in Advanced Chronic Kidney Disease.* N Engl J Med. 2022;387(22):2021-32. — https://pubmed.ncbi.nlm.nih.gov/36326117/ — the randomised result that stopping RAASi in advanced CKD did not improve eGFR — a counterexample that recalls that the hazard layer of the model rests on observational studies.
66. **[C]** Leon SJ, Whitlock R, Rigatto C, et al. *Hyperkalemia-Related Discontinuation of Renin-Angiotensin-Aldosterone System Inhibitors and Clinical Outcomes in CKD.* Am J Kidney Dis. 2022;80(2):164-73. — https://pubmed.ncbi.nlm.nih.gov/35085685/

## 10. Diet, colon, co-medication

67. **[C]** Cupisti A, Kovesdy CP, D'Alessandro C, Kalantar-Zadeh K. *Dietary Approach to Recurrent or Chronic Hyperkalaemia in Patients with Decreased Kidney Function.* Nutrients. 2018;10(3):261. — https://pubmed.ncbi.nlm.nih.gov/29495340/
68. **[C]** Bernier-Jean A, Wong G, Saglimbene V, et al. *Dietary Potassium Intake and All-Cause Mortality in Adults Treated with Hemodialysis.* Clin J Am Soc Nephrol. 2021;16(12):1851-61. — https://pubmed.ncbi.nlm.nih.gov/34853064/ — that the evidence for dietary K restriction is weaker than supposed.
69. **[C]** Neal B, Wu Y, Feng X, et al. (SSaSS) *Effect of Salt Substitution on Cardiovascular Events and Death.* N Engl J Med. 2021;385(12):1067-77. — https://pubmed.ncbi.nlm.nih.gov/34459569/ — a low-sodium salt substitute = a KCl load (the map's `saltsub`).
70. **[S]** Sandle GI, Gaiger E, Tapster S, Goodship TH. *Enhanced rectal potassium secretion in chronic renal insufficiency: evidence for large intestinal potassium adaptation in man.* Clin Sci (Lond). 1986;71(4):393-401. — https://pubmed.ncbi.nlm.nih.gov/3757437/
71. **[C]** Perazella MA. *Drug-induced hyperkalemia: old culprits and new offenders.* Am J Med. 2000;109(4):307-14. — https://pubmed.ncbi.nlm.nih.gov/10996582/
72. **[C]** Antoniou T, Gomes T, Juurlink DN, et al. *Trimethoprim-sulfamethoxazole-induced hyperkalemia in patients receiving inhibitors of the renin-angiotensin system.* Arch Intern Med. 2010;170(12):1045-9. — https://pubmed.ncbi.nlm.nih.gov/20585070/
73. **[C]** Neuen BL, Oshima M, Agarwal R, et al. *Sodium-Glucose Cotransporter 2 Inhibitors and Risk of Hyperkalemia in People With Type 2 Diabetes: A Meta-Analysis.* Circulation. 2022;145(19):1460-70. — https://pubmed.ncbi.nlm.nih.gov/35394821/ — the reduction in hyperkalaemia risk with SGLT2i (the model's `SGLT2I` switch).

## 11. Diabetic ketoacidosis — the same number, the opposite store (the partition patient)

74. **[C]** Kitabchi AE, Umpierrez GE, Miles JM, Fisher JN. *Hyperglycemic crises in adult patients with diabetes.* Diabetes Care. 2009;32(7):1335-43. — https://pubmed.ncbi.nlm.nih.gov/19564476/ — serum K in DKA is normal to raised even with a total body K deficit of 3-5 mmol/kg.
75. **[S]** Adrogué HJ, Lederer ED, Suki WN, Eknoyan G. *Determinants of plasma potassium levels in diabetic ketoacidosis.* Medicine (Baltimore). 1986;65(3):163-72. — https://pubmed.ncbi.nlm.nih.gov/3084904/
76. **[S]** Nicolis GL, Kahn T, Sanchez A, Gabrilove JL. *Glucose-induced hyperkalemia in diabetic subjects.* Arch Intern Med. 1981;141(1):49-53. — https://pubmed.ncbi.nlm.nih.gov/7447584/ — the K shift driven by hyperosmolality (the model's `K_GLU`).

## 12. Clinical guidelines

77. **[C]** Kidney Disease: Improving Global Outcomes (KDIGO) CKD Work Group. *KDIGO 2024 Clinical Practice Guideline for the Evaluation and Management of Chronic Kidney Disease.* Kidney Int. 2024;105(4S):S117-S314. — https://pubmed.ncbi.nlm.nih.gov/38490803/
78. **[C]** Clase CM, Carrero JJ, Ellison DH, et al. *Potassium homeostasis and management of dyskalemia in kidney diseases: conclusions from a KDIGO Controversies Conference.* Kidney Int. 2020;97(1):42-61. — https://pubmed.ncbi.nlm.nih.gov/31706619/
79. **[C]** McDonagh TA, Metra M, Adamo M, et al. *2021 ESC Guidelines for the diagnosis and treatment of acute and chronic heart failure.* Eur Heart J. 2021;42(36):3599-3726. — https://pubmed.ncbi.nlm.nih.gov/34447992/
80. **[C]** Heidenreich PA, Bozkurt B, Aguilar D, et al. *2022 AHA/ACC/HFSA Guideline for the Management of Heart Failure.* Circulation. 2022;145(18):e895-e1032. — https://pubmed.ncbi.nlm.nih.gov/35363499/

## 13. QSP modelling methodology

81. **[S]** Elmokadem A, Riggs MM, Baron KT. *Quantitative Systems Pharmacology and Physiologically-Based Pharmacokinetic Modeling With mrgsolve.* CPT Pharmacometrics Syst Pharmacol. 2019;8(12):883-93. — https://pubmed.ncbi.nlm.nih.gov/31652028/
82. **[S]** Guyton AC, Coleman TG, Granger HJ. *Circulation: overall regulation.* Annu Rev Physiol. 1972;34:13-46. — https://pubmed.ncbi.nlm.nih.gov/4334846/ — the classic source for the fluid compartment volumes (V_ECF 0.20·BW, V_ICF 0.36·BW).
83. **[S]** Hallow KM, Gebremichael Y. *A quantitative systems physiology model of renal function and blood pressure regulation: Model description.* CPT Pharmacometrics Syst Pharmacol. 2017;6(6):383-92. — https://pubmed.ncbi.nlm.nih.gov/28548387/ — a structural reference for renal QSP modelling.
84. **[S]** Hallow KM, Gebremichael Y. *A quantitative systems physiology model of renal function and blood pressure regulation: Application in salt-sensitive hypertension.* CPT Pharmacometrics Syst Pharmacol. 2017;6(6):393-400. — https://pubmed.ncbi.nlm.nih.gov/28556624/

---

## Appendix A — Parameter provenance

| Parameter | Value | Type of source | Basis |
|---|---|---|---|
| `ALPHA` | 0.25 | structure + validation | #4, #5 — chronic buffering 224 mmol/(mmol/L), acute 66 |
| `S0` | 75.18 mmol/day | **fitted** | K 4.20 at an eGFR of 100 (#1, #15) |
| `ADAPT_P` | 0.879 | **fitted** | K 5.00 at an eGFR of 20 (#17, #58) |
| `ADAPT_MX` | 6.36 | **fitted** | K 5.30 at an eGFR of 12 (#17, #25) |
| `KI_MRA` | 0.1029 mg/L | **fitted** | RALES ΔK +0.30 (#34) |
| `N_HCO3` | 0.435 | **fitted** | ΔK -0.30 for HCO3 18→24 (#26, #27) |
| `KC_COL` | 1.00 | structure | upregulated colonic secretion in renal failure (#16, #70) |
| `PHIMAX_P` / `D50_PAT` | 0.42 / 12 g | structure | site of action and mechanism (#44); validation from #39 |
| `PHIMAX_S` / `D50_SZC` | 0.45 / 8 g | structure | site of action and mechanism (#43); validation from #41, #42 |
| `EMAX_INS` | 0.120 | consistent with validation | insulin ΔK -0.6 to -1.0 (#9, #50) |
| `EMAX_B2` | 0.070 | consistent with validation | salbutamol ΔK -0.6 to -1.0 (#9) |
| `K_PH` | 0.190 | structure | mineral acidosis ΔK +0.3 per 0.1 pH (#11, #12) |
| `K_GLU` | 0.0025 | structure | hyperglycaemia-induced K shift (#76) |
| `SK_ALDO` | 1.30 | structure | direct stimulation of aldosterone by K (#32) |
| `BK_HI`,`BK_LO` | 0.642 / 0.916 | **fitted** | the U-shaped mortality risk (#56, #57) |
| `LNHR_RD` | -0.248 | literature value | RAASi HR 0.78 (#62, #63) |
| `SLOPE_RD` | 1.40 mL/min/yr | literature value | (#62, #63) |
| `VTH0`,`H_MID`,`S_CA` | -70, -78, 15 | structure | membrane potential and Na channel availability (#55) |

## Appendix B — What the model does not deal with (explicitly out of scope)

The following are clinically important but were **deliberately not included** in this model, each with its
reference alongside. The model must not be used in these areas.

- **Acute kidney injury (AKI on CKD)** — the commonest cause of genuinely severe hyperkalaemia (#49)
- **Rhabdomyolysis · tumour lysis** — release of K from the ICF store itself; this model treats the ICF only as a passive store (#52)
- **Pseudohyperkalaemia** — haemolysis, leucocytosis, thrombocytosis (#52)
- **Dialysis** — an intermittent, very high clearance removal term (#78)
- **Digoxin toxicity** — direct inhibition of the Na-K-ATPase (#71)
- **Exercise-induced transient hyperkalaemia** — interstitial accumulation in working muscle (#5)
- **Children and pregnancy** — the parameters are for a 70 kg adult

## Appendix C — Citation count

**84 primary references and guidelines** in total (with PubMed links), organised into 13 sections.
