# Bile Acid Diarrhoea / Malabsorption — references

> **Every PMID is a real paper, looked up and verified through the NCBI E-utilities.**
> Link format: `https://pubmed.ncbi.nlm.nih.gov/<PMID>/`
> Each section preamble states **which equation and which parameter** of this QSP model that
> group of references supports. The model's claims are marked so as to distinguish what came from the
> literature from what the model derived for itself (⚙ = derived by the model, 📖 = given by the literature).

---

## Contents

| Section | Subject | Papers |
|---|---|---:|
| A | Quantitative physiology of the enterohepatic circulation — pool size · cycles per day · synthesis rate | 14 |
| B | ASBT / ileal transport — the molecular reality of the `φ` axis | 12 |
| C | The FXR–FGF19–CYP7A1 feedback — the `κ` axis and loop gain | 13 |
| D | Pathophysiology and subtypes of bile acid diarrhoea | 12 |
| E | Colonic secretion · motility — the 3 mM threshold and positive feedback | 11 |
| F | Diagnostic tests — SeHCAT · C4 · FGF19 · faecal bile acids | 14 |
| G | Bile acid sequestrants (class A) | 8 |
| H | FXR agonists (class B) | 9 |
| I | ASBT inhibitors (class C) — the mirror image | 8 |
| J | Ileal resection · the 100 cm rule · steatorrhoea | 7 |
| K | The microbiota — deconjugation and 7α-dehydroxylation (class E) | 8 |
| L | Overlap with IBS-D and symptomatic treatment (class D) | 8 |
| M | Bile acid QSP · PBPK modelling and tools | 10 |
| | **Total** | **134** |

---

## A. Quantitative physiology of the enterohepatic circulation — pool size · cycles per day · synthesis rate

> The model's normal steady-state calibration values come from here: total pool 3 g (6,000 µmol),
> 4–6 cycles per day, duodenal delivery 20–40 mmol/day, hepatic synthesis 0.35 g/day.
> Hofmann 1983 in particular is the first enterohepatic circulation simulation model, and the
> direct ancestor of this repository's model.

1. Hofmann AF, et al. *Description and simulation of a physiological pharmacokinetic model for the metabolism and enterohepatic circulation of bile acids in man.* J Clin Invest 1983. [PMID 6682120](https://pubmed.ncbi.nlm.nih.gov/6682120/)
2. Chiang JY. *Bile acid metabolism and signaling.* Compr Physiol 2013. [PMID 23897684](https://pubmed.ncbi.nlm.nih.gov/23897684/)
3. Chiang JYL. *Bile Acid Metabolism in Liver Pathobiology.* Gene Expr 2018. [PMID 29325602](https://pubmed.ncbi.nlm.nih.gov/29325602/)
4. Chiang JY. *Bile acids: regulation of synthesis.* J Lipid Res 2009. [PMID 19346330](https://pubmed.ncbi.nlm.nih.gov/19346330/)
5. Di Ciaula A, et al. *Bile Acid Physiology.* Ann Hepatol 2017. [PMID 29080336](https://pubmed.ncbi.nlm.nih.gov/29080336/)
6. Roberts MS, et al. *Enterohepatic circulation: physiological, pharmacokinetic and clinical implications.* Clin Pharmacokinet 2002. [PMID 12162761](https://pubmed.ncbi.nlm.nih.gov/12162761/)
7. Vantrappen G, et al. *A new method for the measurement of bile acid turnover and pool size by a double label, single intubation technique.* J Lipid Res 1981. [PMID 7017051](https://pubmed.ncbi.nlm.nih.gov/7017051/)
8. Hulzebos CV, et al. *Measurement of parameters of cholic acid kinetics in plasma using a microscale stable isotope dilution technique.* J Lipid Res 2001. [PMID 11714862](https://pubmed.ncbi.nlm.nih.gov/11714862/)
9. Stellaard F, et al. *Determination of deoxycholic acid pool size and input rate using [24-13C]deoxycholic acid.* J Lipid Res 1986. [PMID 3559388](https://pubmed.ncbi.nlm.nih.gov/3559388/)
10. Schwartz CC, et al. *Cholesterol kinetics in subjects with bile fistula. Positive relationship between size of the bile acid precursor pool and bile acid synthetic rate.* J Clin Invest 1993. [PMID 8450070](https://pubmed.ncbi.nlm.nih.gov/8450070/)
11. Balistreri WF, et al. *Validation of use of 11,12-2H-labeled chenodeoxycholic acid in isotope dilution measurements of bile acid kinetics.* Pediatr Res 1975. [PMID 1103070](https://pubmed.ncbi.nlm.nih.gov/1103070/)
12. Javitt NB. *Bile acid synthesis from cholesterol: regulatory and auxiliary pathways.* FASEB J 1994. [PMID 8001744](https://pubmed.ncbi.nlm.nih.gov/8001744/)
13. Schwarz M, et al. *Two 7 alpha-hydroxylase enzymes in bile acid biosynthesis.* Curr Opin Lipidol 1998. [PMID 9559267](https://pubmed.ncbi.nlm.nih.gov/9559267/)
14. Sonne DP, et al. *Postprandial gallbladder emptying in patients with type 2 diabetes.* Eur J Endocrinol 2014. [PMID 24986531](https://pubmed.ncbi.nlm.nih.gov/24986531/)

---

## B. ASBT / ileal transport — the molecular reality of the `φ` axis

> The model's first lesion axis `φ` is the **surviving fraction of ileal ASBT absorptive capacity**.
> The SLC10A2 mutation of Oelkers 1997 is the human experiment for `φ → 0`, and
> the Slc10a2 knockout of Dawson 2003 is its animal counterpart. Hruz 2006 showed in humans
> that ASBT really is reduced in ileal disease.

15. Oelkers P, et al. *Primary bile acid malabsorption caused by mutations in the ileal sodium-dependent bile acid transporter gene (SLC10A2).* J Clin Invest 1997. [PMID 9109432](https://pubmed.ncbi.nlm.nih.gov/9109432/)
16. Dawson PA, et al. *Targeted deletion of the ileal bile acid transporter eliminates enterohepatic cycling of bile acids in mice.* J Biol Chem 2003. [PMID 12819193](https://pubmed.ncbi.nlm.nih.gov/12819193/)
17. Dawson PA, et al. *Bile acid transporters.* J Lipid Res 2009. [PMID 19498215](https://pubmed.ncbi.nlm.nih.gov/19498215/)
18. Hruz P, et al. *Adaptive regulation of the ileal apical sodium dependent bile acid transporter (ASBT) in patients with obstructive cholestasis.* Gut 2006. [PMID 16150853](https://pubmed.ncbi.nlm.nih.gov/16150853/)
19. Balesaria S, et al. *Exploring possible mechanisms for primary bile acid malabsorption: evidence for different regulation of ileal bile acid transporter transcripts.* Eur J Gastroenterol Hepatol 2008. [PMID 18403943](https://pubmed.ncbi.nlm.nih.gov/18403943/)
20. Renner O, et al. *Mutation screening of apical sodium-dependent bile acid transporter (SLC10A2).* Hum Genet 2009. [PMID 19184108](https://pubmed.ncbi.nlm.nih.gov/19184108/)
21. Dawson PA, et al. *The heteromeric organic solute transporter alpha-beta, Ostalpha-Ostbeta, is an ileal basolateral bile acid transporter.* J Biol Chem 2005. [PMID 15563450](https://pubmed.ncbi.nlm.nih.gov/15563450/)
22. Rao A, et al. *The organic solute transporter alpha-beta, Ostalpha-Ostbeta, is essential for intestinal bile acid transport and homeostasis.* PNAS 2008. [PMID 18292224](https://pubmed.ncbi.nlm.nih.gov/18292224/)
23. Ballatori N, et al. *OSTalpha-OSTbeta: a major basolateral bile acid and steroid transporter.* Hepatology 2005. [PMID 16317684](https://pubmed.ncbi.nlm.nih.gov/16317684/)
24. Sultan M, et al. *Organic solute transporter-β (SLC51B) deficiency in two brothers with congenital diarrhea.* Hepatology 2018. [PMID 28898457](https://pubmed.ncbi.nlm.nih.gov/28898457/)
25. Wang L, et al. *Mechanism of Asbt (Slc10a2)-related bile acid malabsorption in diarrhea after pelvic radiation.* Int J Radiat Biol 2020. [PMID 31900034](https://pubmed.ncbi.nlm.nih.gov/31900034/)
26. Chothe PP, et al. *Human bile acid transporter ASBT (SLC10A2) forms functional non-covalent homodimers.* Biochim Biophys Acta 2018. [PMID 29198943](https://pubmed.ncbi.nlm.nih.gov/29198943/)

---

## C. The FXR–FGF19–CYP7A1 feedback — the `κ` axis and loop gain

> This model's central claim — **the sensor measures the absorbed flux while the symptoms are made by the
> flux that overflows** — comes from this group of references. Inagaki 2005 established that ileal
> FGF15/19 is the endocrine signal that suppresses CYP7A1, and Lundåsen 2006
> showed a marked diurnal variation of FGF19 in humans and its control of bile acid synthesis.
> Lu 2000 / Kerr 2002 showed that the SHP arm is **auxiliary** (feedback survives its removal),
> and that is the basis for giving the FGF19 arm the greater weight (Rf dominates) in the model.

27. Inagaki T, et al. *Fibroblast growth factor 15 functions as an enterohepatic signal to regulate bile acid homeostasis.* Cell Metab 2005. [PMID 16213224](https://pubmed.ncbi.nlm.nih.gov/16213224/)
28. Lundåsen T, et al. *Circulating intestinal fibroblast growth factor 19 has a pronounced diurnal variation and modulates hepatic bile acid synthesis in man.* J Intern Med 2006. [PMID 17116003](https://pubmed.ncbi.nlm.nih.gov/17116003/)
29. Lu TT, et al. *Molecular basis for feedback regulation of bile acid synthesis by nuclear receptors.* Mol Cell 2000. [PMID 11030331](https://pubmed.ncbi.nlm.nih.gov/11030331/)
30. Kerr TA, et al. *Loss of nuclear receptor SHP impairs but does not eliminate negative feedback regulation of bile acid synthesis.* Dev Cell 2002. [PMID 12062084](https://pubmed.ncbi.nlm.nih.gov/12062084/)
31. Makishima M, et al. *Identification of a nuclear receptor for bile acids.* Science 1999. [PMID 10334992](https://pubmed.ncbi.nlm.nih.gov/10334992/)
32. Fu T, et al. *FXR Primes the Liver for Intestinal FGF15 Signaling by Transient Induction of β-Klotho.* Mol Endocrinol 2016. [PMID 26505219](https://pubmed.ncbi.nlm.nih.gov/26505219/)
33. Byun S, et al. *Postprandial FGF19-induced phosphorylation by Src is critical for FXR function in bile acid homeostasis.* Nat Commun 2018. [PMID 29968724](https://pubmed.ncbi.nlm.nih.gov/29968724/)
34. Wang C, et al. *Hepatocyte FRS2α is essential for the endocrine fibroblast growth factor to limit the amplitude of bile acid production.* Curr Mol Med 2014. [PMID 25056539](https://pubmed.ncbi.nlm.nih.gov/25056539/)
35. Bouju A, et al. *A primer on the pleiotropic endocrine fibroblast growth factor FGF19/FGF15.* Differentiation 2024. [PMID 39500656](https://pubmed.ncbi.nlm.nih.gov/39500656/)
36. Ferrell JM, Chiang JYL. *Circadian rhythms in liver metabolism and disease.* Acta Pharm Sin B 2015. [PMID 26579436](https://pubmed.ncbi.nlm.nih.gov/26579436/)
37. Al-Khaifi A, et al. *Asynchronous rhythms of circulating conjugated and unconjugated bile acids.* J Intern Med 2018. [PMID 29964306](https://pubmed.ncbi.nlm.nih.gov/29964306/)
38. Haeusler RA, et al. *Impaired generation of 12-hydroxylated bile acids links hepatic insulin signaling with dyslipidemia.* Cell Metab 2012. [PMID 22197325](https://pubmed.ncbi.nlm.nih.gov/22197325/)
39. Patankar JV, et al. *Genetic ablation of Cyp8b1 preserves host metabolic function.* Am J Physiol Endocrinol Metab 2018. [PMID 29066462](https://pubmed.ncbi.nlm.nih.gov/29066462/)

---

## D. Pathophysiology and subtypes of bile acid diarrhoea

> Type 1 (ileal disease/resection, `φ↓`) · Type 2 (idiopathic/primary, `κ↓`) ·
> Type 3 (cholecystectomy · SIBO · coeliac disease and the like) · Type 4 (drug-induced).
> The Walters 2009 line of work, showing that **fasting FGF19 is low** in idiopathic BAD, is
> the direct basis for introducing the `κ` axis (see Pattni 2012 and
> Borup 2015 in section F).

40. Camilleri M. *Bile Acid Diarrhea: Prevalence, Pathogenesis, and Therapy* / *The Role of Bile Acids in Chronic Diarrhea.* Am J Gastroenterol 2020. [PMID 32558690](https://pubmed.ncbi.nlm.nih.gov/32558690/)
41. Camilleri M, Vijayvargiya P. *Bile Acid Diarrhea in Adults and Adolescents.* Neurogastroenterol Motil 2022. [PMID 34751982](https://pubmed.ncbi.nlm.nih.gov/34751982/)
42. Walters JRF. *Managing bile acid diarrhea: aspects of contention.* Expert Rev Gastroenterol Hepatol 2024. [PMID 39264409](https://pubmed.ncbi.nlm.nih.gov/39264409/)
43. Barbara G, et al. *Bile acid diarrhea in patients with chronic diarrhea. Current appraisal and recommendations for clinical practice.* Dig Liver Dis 2025. [PMID 39827025](https://pubmed.ncbi.nlm.nih.gov/39827025/)
44. Di Ciaula A, et al. *Advances in the pathophysiology, diagnosis and management of chronic diarrhoea from bile acid malabsorption.* Eur J Intern Med 2024. [PMID 39069430](https://pubmed.ncbi.nlm.nih.gov/39069430/)
45. Barkun AN, et al. *Bile acid malabsorption in chronic diarrhea: pathophysiology and treatment.* Can J Gastroenterol 2013. [PMID 24199211](https://pubmed.ncbi.nlm.nih.gov/24199211/)
46. Gracie DJ, et al. *Prevalence of, and predictors of, bile acid malabsorption in outpatients with chronic diarrhea.* Neurogastroenterol Motil 2012. [PMID 22765392](https://pubmed.ncbi.nlm.nih.gov/22765392/)
47. Ruiz-Campos L, et al. *Systematic review with meta-analysis: the prevalence of bile acid malabsorption and response to colestyramine.* Aliment Pharmacol Ther 2019. [PMID 30585336](https://pubmed.ncbi.nlm.nih.gov/30585336/)
48. Farrugia A, et al. *Rates of Bile Acid Diarrhoea After Cholecystectomy: A Multicentre Audit.* World J Surg 2021. [PMID 33982189](https://pubmed.ncbi.nlm.nih.gov/33982189/)
49. Farrugia A, et al. *Bile acid diarrhoea and metabolic changes after cholecystectomy: a prospective case-control study.* BMC Gastroenterol 2024. [PMID 39174936](https://pubmed.ncbi.nlm.nih.gov/39174936/)
50. Barrera F, et al. *Effect of cholecystectomy on bile acid synthesis and circulating levels of fibroblast growth factor 19.* Ann Hepatol 2016. [PMID 26256900](https://pubmed.ncbi.nlm.nih.gov/26256900/)
51. Bouchoucha M, et al. *Metformin and digestive disorders.* Diabetes Metab 2011. [PMID 21236717](https://pubmed.ncbi.nlm.nih.gov/21236717/)

---

## E. Colonic secretion · motility — the 3 mM threshold and positive feedback

> The model's `EC50_G = 3 mM CDCA-equivalent` comes directly from the human colonic perfusion
> experiment of **Mekjian 1971**: above 3 mM, the dihydroxy bile acids (CDCA·DCA) switch the colon
> to net secretion. The trihydroxy acid (CA) is far weaker → hence the model's potency
> weights (CA 0.15 / CDCA 1.00 / DCA 1.15 / LCA 0.30).
> Farack 1984 showed that loperamide inhibits DCA-induced secretion itself,
> and that is the basis on which the model classifies loperamide not as "symptomatic" but as a
> **loop-gain reducer**.

52. Mekjian HS, Phillips SF, Hofmann AF. *Colonic secretion of water and electrolytes induced by bile acids: perfusion studies in man.* J Clin Invest 1971. [PMID 4938344](https://pubmed.ncbi.nlm.nih.gov/4938344/)
53. Rao AS, et al. *Chenodeoxycholate in females with irritable bowel syndrome-constipation: a pharmacodynamic and pharmacogenetic analysis.* Gastroenterology 2010. [PMID 20691689](https://pubmed.ncbi.nlm.nih.gov/20691689/)
54. Odunsi-Shiyanbade ST, et al. *Effects of chenodeoxycholate and a bile acid sequestrant, colesevelam, on intestinal transit and bowel function.* Clin Gastroenterol Hepatol 2010. [PMID 19879973](https://pubmed.ncbi.nlm.nih.gov/19879973/)
55. Alemi F, et al. *The receptor TGR5 mediates the prokinetic actions of intestinal bile acids and is required for normal defecation.* Gastroenterology 2013. [PMID 23041323](https://pubmed.ncbi.nlm.nih.gov/23041323/)
56. Ward JB, et al. *The bile acid receptor, TGR5, regulates basal and cholinergic-induced secretory responses in rat colon.* Neurogastroenterol Motil 2013. [PMID 23634890](https://pubmed.ncbi.nlm.nih.gov/23634890/)
57. Poole DP, et al. *Expression and function of the bile acid receptor GpBAR1 (TGR5) in the murine enteric nervous system.* Neurogastroenterol Motil 2010. [PMID 20236244](https://pubmed.ncbi.nlm.nih.gov/20236244/)
58. Bunnett NW. *Neuro-humoral signalling by bile acids and the TGR5 receptor in the gastrointestinal tract.* J Physiol 2014. [PMID 24614746](https://pubmed.ncbi.nlm.nih.gov/24614746/)
59. Yu Y, et al. *Deoxycholic acid activates colonic afferent nerves via 5-HT3 receptor-dependent and -independent mechanisms.* Am J Physiol Gastrointest Liver Physiol 2019. [PMID 31216174](https://pubmed.ncbi.nlm.nih.gov/31216174/)
60. Camilleri M, et al. *Bile acid detergency: permeability, inflammation, and effects of sulfation.* Am J Physiol Gastrointest Liver Physiol 2022. [PMID 35258349](https://pubmed.ncbi.nlm.nih.gov/35258349/)
61. Zeng H, et al. *Deoxycholic Acid Modulates Cell-Junction Gene Expression and Increases Intestinal Barrier Dysfunction.* Molecules 2022. [PMID 35163990](https://pubmed.ncbi.nlm.nih.gov/35163990/)
62. Farack UM, Loeschke K. *Inhibition by loperamide of deoxycholic acid induced intestinal secretion.* Naunyn Schmiedebergs Arch Pharmacol 1984. [PMID 6328335](https://pubmed.ncbi.nlm.nih.gov/6328335/)

---

## F. Diagnostic tests — SeHCAT · C4 · FGF19 · faecal bile acids

> **The model's most testable prediction either conflicts with or agrees with this section.**
> ⚙ The model claims that SeHCAT reads `φ` alone (and reads it compressed into a very narrow range), while
> C4 and faecal bile acids read `S = f(φ, κ)`. Hence the two tests must agree only
> imperfectly (the Valentin 2016 meta-analysis), and patients must exist with a normal SeHCAT but
> a raised C4 (Vijayvargiya 2017·2020, Borup 2024).

63. Riemsma R, et al. *SeHCAT [tauroselcholic (selenium-75) acid] for the investigation of bile acid malabsorption: a systematic review and cost-effectiveness analysis.* Health Technol Assess 2013. [PMID 24351663](https://pubmed.ncbi.nlm.nih.gov/24351663/)
64. Wedlake L, et al. *Systematic review: the prevalence of idiopathic bile acid malabsorption as diagnosed by SeHCAT scanning in patients with diarrhoea-predominant IBS.* Aliment Pharmacol Ther 2009. [PMID 19570102](https://pubmed.ncbi.nlm.nih.gov/19570102/)
65. van Tilburg AJ, et al. *The selenium-75-homocholic acid taurine test reevaluated: combined measurement of fecal selenium-75 activity and 3 alpha-hydroxy bile acids.* J Nucl Med 1991. [PMID 2045936](https://pubmed.ncbi.nlm.nih.gov/2045936/)
66. Pattni SS, et al. *Fibroblast Growth Factor 19 and 7α-Hydroxy-4-Cholesten-3-one in the Diagnosis of Patients With Possible Bile Acid Diarrhea.* Clin Transl Gastroenterol 2012. [PMID 23238290](https://pubmed.ncbi.nlm.nih.gov/23238290/)
67. Vijayvargiya P, et al. *Performance characteristics of serum C4 and FGF19 measurements to exclude the diagnosis of bile acid diarrhoea.* Aliment Pharmacol Ther 2017. [PMID 28691284](https://pubmed.ncbi.nlm.nih.gov/28691284/)
68. Vijayvargiya P, et al. *Fecal Bile Acid Testing in Assessing Patients With Chronic Unexplained Diarrhea.* Am J Gastroenterol 2020. [PMID 32618660](https://pubmed.ncbi.nlm.nih.gov/32618660/)
69. Valentin N, et al. *Biomarkers for bile acid diarrhoea in functional bowel disorder with diarrhoea: a systematic review and meta-analysis.* Gut 2016. [PMID 26347530](https://pubmed.ncbi.nlm.nih.gov/26347530/)
70. Borup C, et al. *Diagnosis of bile acid diarrhoea by fasting and postprandial measurements of fibroblast growth factor 19.* Eur J Gastroenterol Hepatol 2015. [PMID 26426834](https://pubmed.ncbi.nlm.nih.gov/26426834/)
71. Borup C, et al. *Prospective comparison of diagnostic tests for bile acid diarrhoea.* Aliment Pharmacol Ther 2024. [PMID 37794830](https://pubmed.ncbi.nlm.nih.gov/37794830/)
72. Battat R, et al. *Serum Concentrations of 7α-hydroxy-4-cholesten-3-one Are Associated With Bile Acid Diarrhea in Patients With Crohn's Disease.* Clin Gastroenterol Hepatol 2019. [PMID 30448597](https://pubmed.ncbi.nlm.nih.gov/30448597/)
73. Lyutakov I, et al. *Methods for diagnosing bile acid malabsorption: a systematic review.* BMC Gastroenterol 2019. [PMID 31726982](https://pubmed.ncbi.nlm.nih.gov/31726982/)
74. Dilmaghani S, et al. *Simplifying Diagnosis of Bile Acid Diarrhea With Clinical and Biochemical Measurements on Blood and Single Stool Sample.* Clin Gastroenterol Hepatol 2025. [PMID 40378992](https://pubmed.ncbi.nlm.nih.gov/40378992/)
75. Peleman C, et al. *Colonic Transit and Bile Acid Synthesis or Excretion in Patients With Irritable Bowel Syndrome-Diarrhea Without Bile Acid Malabsorption.* Clin Gastroenterol Hepatol 2017. [PMID 27856362](https://pubmed.ncbi.nlm.nih.gov/27856362/)
76. Lupianez-Merly C, Camilleri M. *Recent developments in diagnosing bile acid diarrhea.* Expert Rev Gastroenterol Hepatol 2023. [PMID 38086533](https://pubmed.ncbi.nlm.nih.gov/38086533/)

---

## G. Bile acid sequestrants (class A — they bind the luminal load)

> ⚙ The model predicts that a sequestrant **raises the very driving force it is treating**:
> binding luminal bile acids removes the sensor ligand too, so FGF19 falls and
> CYP7A1/C4 rise. This "structural escape" is why the dose-response of a sequestrant is
> shallow, and why combination with an FXR agonist is synergistic.
> 📖 That a sequestrant really does increase bile acid synthesis is a long-established fact in the
> lipid literature (it is the LDL-lowering mechanism itself).

77. Borup C, et al. *Efficacy and safety of colesevelam for the treatment of bile acid diarrhoea: a double-blind, randomised, placebo-controlled, phase 4 clinical trial.* Lancet Gastroenterol Hepatol 2023. [PMID 36758570](https://pubmed.ncbi.nlm.nih.gov/36758570/)
78. Camilleri M, et al. *Effect of colesevelam on faecal bile acids and bowel functions in diarrhoea-predominant irritable bowel syndrome.* Aliment Pharmacol Ther 2015. [PMID 25594801](https://pubmed.ncbi.nlm.nih.gov/25594801/)
79. Beigel F, et al. *Colesevelam for the treatment of bile acid malabsorption-associated diarrhea in patients with Crohn's disease.* J Crohns Colitis 2014. [PMID 24953836](https://pubmed.ncbi.nlm.nih.gov/24953836/)
80. Soares GAR, et al. *Efficacy of Bile Acid Sequestrants in the Treatment of Bile Acid Diarrhea: A Meta-Analysis of Randomized Controlled Trials.* J Clin Pharmacol 2025. [PMID 39428959](https://pubmed.ncbi.nlm.nih.gov/39428959/)
81. Wilcox C, et al. *Systematic review: the management of chronic diarrhoea due to bile acid malabsorption.* Aliment Pharmacol Ther 2014. [PMID 24602022](https://pubmed.ncbi.nlm.nih.gov/24602022/)
82. Kårhus ML, et al. *Safety and efficacy of liraglutide versus colesevelam for the treatment of bile acid diarrhoea: a randomised, double-blind, active-comparator, non-inferiority clinical trial.* Lancet Gastroenterol Hepatol 2022. [PMID 35868334](https://pubmed.ncbi.nlm.nih.gov/35868334/)
83. Ellegaard AM, et al. *Liraglutide and Colesevelam Change Serum and Fecal Bile Acid Levels in a Randomized Trial.* Clin Transl Gastroenterol 2024. [PMID 39602188](https://pubmed.ncbi.nlm.nih.gov/39602188/)
84. Scaldaferri F, et al. *Use and indications of cholestyramine and bile acid sequestrants.* Intern Emerg Med 2013. [PMID 21739227](https://pubmed.ncbi.nlm.nih.gov/21739227/)

---

## H. FXR agonists (class B — they supply the missing signal)

> ⚙ The model's classification: an FXR agonist is **additive** in receptor occupancy, so it works
> independently of `φ`, and it is the only class that lowers `S` (= the colonic load) itself.
> ⚙ Prediction: an FXR agonist raises the SeHCAT value without fixing the ileum
> (because the pool turnover `S/P` falls). 📖 Walters 2015 (OBADIAH) and
> Camilleri 2020 (tropifexor) showed a sharp fall in C4, a rise in FGF19 and slowed colonic transit
> in humans. 📖 The price is a rise in LDL-C (Siddiqui 2020, Pencek 2016).

85. Walters JR, et al. *The response of patients with bile acid diarrhoea to the farnesoid X receptor agonist obeticholic acid.* Aliment Pharmacol Ther 2015. [PMID 25329562](https://pubmed.ncbi.nlm.nih.gov/25329562/)
86. Camilleri M, et al. *Randomised clinical trial: significant biochemical and colonic transit effects of the farnesoid X receptor agonist tropifexor in patients with primary bile acid diarrhoea.* Aliment Pharmacol Ther 2020. [PMID 32702169](https://pubmed.ncbi.nlm.nih.gov/32702169/)
87. Pellicciari R, et al. *6alpha-ethyl-chenodeoxycholic acid (6-ECDCA), a potent and selective FXR agonist.* J Med Chem 2002. [PMID 12166927](https://pubmed.ncbi.nlm.nih.gov/12166927/)
88. Jiang L, et al. *Structural basis of tropifexor as a potent and selective agonist of farnesoid X receptor.* Biochem Biophys Res Commun 2021. [PMID 33121679](https://pubmed.ncbi.nlm.nih.gov/33121679/)
89. Gege C, et al. *Nonsteroidal FXR Ligands: Current Status and Clinical Applications.* Handb Exp Pharmacol 2019. [PMID 31197565](https://pubmed.ncbi.nlm.nih.gov/31197565/)
90. Fiorucci S, et al. *Obeticholic Acid: An Update of Its Pharmacological Activities in Liver Disorders.* Handb Exp Pharmacol 2019. [PMID 31201552](https://pubmed.ncbi.nlm.nih.gov/31201552/)
91. Siddiqui MS, et al. *Impact of obeticholic acid on the lipoprotein profile in patients with non-alcoholic steatohepatitis.* J Hepatol 2020. [PMID 31634532](https://pubmed.ncbi.nlm.nih.gov/31634532/)
92. Pencek R, et al. *Effects of obeticholic acid on lipoprotein metabolism in healthy volunteers.* Diabetes Obes Metab 2016. [PMID 27109453](https://pubmed.ncbi.nlm.nih.gov/27109453/)
93. Kremoser C. *FXR agonists for NASH: How are they different and what difference do they make?* J Hepatol 2021. [PMID 33985820](https://pubmed.ncbi.nlm.nih.gov/33985820/)

---

## I. ASBT inhibitors (class C — they lower `φ` on purpose) · the mirror image

> ⚙ The model's **falsifiability test**: an ASBT inhibitor ought to produce the full BAM
> phenotype in a normal intestine, and that ought to be why elobixibat is a treatment for chronic
> constipation. 📖 The Simrén 2011 and Nakajima/Miner line of trials are exactly that.
> 📖 Carreño 2025 is a pharmacokinetic-pharmacodynamic analysis **predicting IBAT-inhibitor-related
> diarrhoea from the C4 concentration**, which uses in clinical trial data exactly the same logic as
> this model's axis that "C4 reads `S`".

94. Simrén M, et al. *Randomised clinical trial: the ileal bile acid transporter inhibitor A3309 vs. placebo in patients with chronic idiopathic constipation.* Aliment Pharmacol Ther 2011. [PMID 21545606](https://pubmed.ncbi.nlm.nih.gov/21545606/)
95. Miner PB Jr. *Elobixibat, the first-in-class Ileal Bile Acid Transporter inhibitor, for the treatment of Chronic Idiopathic Constipation.* Expert Opin Pharmacother 2018. [PMID 30129377](https://pubmed.ncbi.nlm.nih.gov/30129377/)
96. Taniguchi S, et al. *Elobixibat, an ileal bile acid transporter inhibitor, induces giant migrating contractions during natural defecation.* Neurogastroenterol Motil 2018. [PMID 30129138](https://pubmed.ncbi.nlm.nih.gov/30129138/)
97. Manabe N, et al. *Elobixibat Improves Stool/Gas Distribution and Fecal Bile Acids in Older Adults With Chronic Constipation.* JGH Open 2025. [PMID 40832007](https://pubmed.ncbi.nlm.nih.gov/40832007/)
98. Carreño F, et al. *Analysis of C4 Concentrations to Predict Impact of Patient-Reported Diarrhea Associated With the Ileal Bile Acid Transporter Inhibitor...* CPT Pharmacometrics Syst Pharmacol 2025. [PMID 39945351](https://pubmed.ncbi.nlm.nih.gov/39945351/)
99. Al-Dury S, Marschall HU. *Ileal Bile Acid Transporter Inhibition for the Treatment of Chronic Constipation, Cholestatic Pruritus, and NASH.* Front Pharmacol 2018. [PMID 30186169](https://pubmed.ncbi.nlm.nih.gov/30186169/)
100. Shirley M. *Maralixibat: First Approval.* Drugs 2022. [PMID 34813049](https://pubmed.ncbi.nlm.nih.gov/34813049/)
101. Peverelle M, et al. *Review Article: Ileal Bile Acid Transport (IBAT) Inhibitors as an Emerging Treatment for Cholestatic Liver Disease.* Aliment Pharmacol Ther 2026. [PMID 41953994](https://pubmed.ncbi.nlm.nih.gov/41953994/)

---

## J. Ileal resection · the 100 cm rule · steatorrhoea

> ⚙ The model **derives the 100 cm rule rather than fitting it to data**:
> the distal weighting of ASBT density (decay length 60 cm) + the CYP7A1 compensation ceiling (6-fold) +
> the critical micellar concentration (1.5 mM) — these three alone give the transition point from
> bile acid diarrhoea to steatorrhoea. 📖 Hofmann & Poley 1972 / Poley & Hofmann 1976 are the clinical
> primary sources, and Skouras 2019 showed the correlation between resection length and BAM severity.

102. Hofmann AF, Poley JR. *Role of bile acid malabsorption in pathogenesis of diarrhea and steatorrhea in patients with ileal resection.* Gastroenterology 1972. [PMID 5029077](https://pubmed.ncbi.nlm.nih.gov/5029077/)
103. Poley JR, Hofmann AF. *Role of fat maldigestion in pathogenesis of steatorrhea in ileal resection.* Gastroenterology 1976. [PMID 6360](https://pubmed.ncbi.nlm.nih.gov/6360/)
104. Skouras T, et al. *Brief report: length of ileal resection correlates with severity of bile acid malabsorption in Crohn's disease.* Int J Colorectal Dis 2019. [PMID 30116880](https://pubmed.ncbi.nlm.nih.gov/30116880/)
105. Akerlund JE, et al. *Hepatic metabolism of cholesterol in Crohn's disease. Effect of partial resection of ileum.* Gastroenterology 1991. [PMID 2001802](https://pubmed.ncbi.nlm.nih.gov/2001802/)
106. Färkkilä MA, et al. *Plasma lathosterol as a screening test for bile acid malabsorption due to ileal resection.* Clin Sci (Lond) 1996. [PMID 8777839](https://pubmed.ncbi.nlm.nih.gov/8777839/)
107. Siener R, et al. *Intestinal Oxalate Absorption, Enteric Hyperoxaluria, and Risk of Urinary Stone Formation.* Nutrients 2024. [PMID 38257157](https://pubmed.ncbi.nlm.nih.gov/38257157/)
108. Larsen HM, et al. *Chronic loose stools following right-sided hemicolectomy for colon cancer and the association with bile acid malabsorption.* Colorectal Dis 2023. [PMID 36347822](https://pubmed.ncbi.nlm.nih.gov/36347822/)

---

## K. The microbiota — deconjugation and 7α-dehydroxylation (class E — potency, not load)

> ⚙ The model's second account of a discrepancy: the microbiota (almost) does not change the
> **amount** of bile acid reaching the colon but changes its **secretory potency**. Because the
> 7α-dehydroxylation that turns CA (w 0.15) into DCA (w 1.15) is the rate-limiting step, symptoms can
> differ threefold at one and the same faecal bile acid excretion.

109. Ridlon JM, et al. *Bile salt biotransformations by human intestinal bacteria.* J Lipid Res 2006. [PMID 16299351](https://pubmed.ncbi.nlm.nih.gov/16299351/)
110. Wells JE, Hylemon PB. *Identification and characterization of a bile acid 7alpha-dehydroxylation operon in Clostridium sp. strain TO-931.* Appl Environ Microbiol 2000. [PMID 10698778](https://pubmed.ncbi.nlm.nih.gov/10698778/)
111. Ridlon JM, et al. *Identification and characterization of two bile acid coenzyme A transferases from Clostridium scindens.* J Lipid Res 2012. [PMID 22021638](https://pubmed.ncbi.nlm.nih.gov/22021638/)
112. Ridlon JM, et al. *Consequences of bile salt biotransformations by intestinal bacteria.* Gut Microbes 2016. [PMID 26939849](https://pubmed.ncbi.nlm.nih.gov/26939849/)
113. Ridlon JM, et al. *The human gut sterolbiome: bile acid-microbiome endocrine aspects and therapeutics.* Acta Pharm Sin B 2015. [PMID 26579434](https://pubmed.ncbi.nlm.nih.gov/26579434/)
114. Guzior DV, et al. *Bile salt hydrolase acyltransferase activity expands bile acid diversity.* Nature 2024. [PMID 38326608](https://pubmed.ncbi.nlm.nih.gov/38326608/)
115. Sinha SR, et al. *Dysbiosis-Induced Secondary Bile Acid Deficiency Promotes Intestinal Inflammation.* Cell Host Microbe 2020. [PMID 32101703](https://pubmed.ncbi.nlm.nih.gov/32101703/)
116. Camilleri M, et al. *Comparison of biochemical, microbial and mucosal mRNA expression in bile acid diarrhoea and irritable bowel syndrome with diarrhoea.* Gut 2023. [PMID 35580964](https://pubmed.ncbi.nlm.nih.gov/35580964/)

---

## L. Overlap with IBS-D and symptomatic treatment (class D — they lower loop gain)

> ⚙ In the model, loperamide and ondansetron are not merely symptomatic therapy:
> slowing transit lengthens ileal contact time, so `f` rises and the spillover itself
> falls. The model therefore predicts that under loperamide **faecal bile acid excretion too**
> (not only water) will fall.

117. Camilleri M. *Irritable bowel syndrome: treatment based on pathophysiology and biomarkers.* Gut 2023. [PMID 36307180](https://pubmed.ncbi.nlm.nih.gov/36307180/)
118. Garsed K, et al. *A randomised trial of ondansetron for the treatment of irritable bowel syndrome with diarrhoea.* Gut 2014. [PMID 24334242](https://pubmed.ncbi.nlm.nih.gov/24334242/)
119. Gunn D, et al. *Randomised, placebo-controlled trial and meta-analysis show benefit of ondansetron for irritable bowel syndrome with diarrhoea (TRITON).* Aliment Pharmacol Ther 2023. [PMID 36866724](https://pubmed.ncbi.nlm.nih.gov/36866724/)
120. Sun WM, et al. *Effects of loperamide oxide on gastrointestinal transit time and anorectal function in patients with chronic diarrhoea.* Scand J Gastroenterol 1997. [PMID 9018764](https://pubmed.ncbi.nlm.nih.gov/9018764/)
121. Mainguet P, Fiasse R. *Double-blind placebo-controlled study of loperamide (Imodium) in chronic diarrhoea caused by ileocolic disease or resection.* Gut 1977. [PMID 326642](https://pubmed.ncbi.nlm.nih.gov/326642/)
122. Vijayvargiya P, et al. *Safety and Efficacy of Eluxadoline in Patients with Irritable Bowel Syndrome-Diarrhea With or Without Bile Acid Diarrhea.* Dig Dis Sci 2022. [PMID 35122592](https://pubmed.ncbi.nlm.nih.gov/35122592/)
123. Shin A, et al. *Bowel functions, fecal unconjugated primary and secondary bile acids, and colonic transit in patients with irritable bowel syndrome.* Clin Gastroenterol Hepatol 2013. [PMID 23639599](https://pubmed.ncbi.nlm.nih.gov/23639599/)
124. Vijayvargiya P, et al. *Bile and fat excretion are biomarkers of clinically significant diarrhoea and constipation in irritable bowel syndrome.* Aliment Pharmacol Ther 2019. [PMID 30740753](https://pubmed.ncbi.nlm.nih.gov/30740753/)

---

## M. Bile acid QSP · PBPK modelling and tools

> The methodological ancestors of this model. Hofmann 1983 (section A-1) is the prototype,
> Voronova 2020 the modern PBPK reconstruction, Guiastrennec 2018 an enterohepatic circulation model of
> human plasma bile acids, and Meessen 2020 a model of the individual postprandial response.

125. Voronova V, et al. *A Physiology-Based Model of Bile Acid Distribution and Metabolism Under Healthy and Pathologic Conditions.* Cell Mol Gastroenterol Hepatol 2020. [PMID 32112828](https://pubmed.ncbi.nlm.nih.gov/32112828/)
126. Guiastrennec B, et al. *Model-Based Prediction of Plasma Concentration and Enterohepatic Circulation of Total Bile Acids in Humans.* CPT Pharmacometrics Syst Pharmacol 2018. [PMID 30070437](https://pubmed.ncbi.nlm.nih.gov/30070437/)
127. Meessen ECE, et al. *Model-based data analysis of individual human postprandial plasma bile acid responses.* Physiol Rep 2020. [PMID 32170845](https://pubmed.ncbi.nlm.nih.gov/32170845/)
128. Mayo AK, et al. *Quantitative Systems Toxicology Model Predicts Obeticholic Acid-Associated Liver Injury.* Clin Pharmacol Ther 2026. [PMID 42332345](https://pubmed.ncbi.nlm.nih.gov/42332345/)
129. Generaux G, et al. *Quantitative systems toxicology (QST) reproduces species differences in PF-04895162 liver safety due to combined mitochondrial and bile acid toxicity.* Pharmacol Res Perspect 2019. [PMID 31624633](https://pubmed.ncbi.nlm.nih.gov/31624633/)
130. Nigam SK. *The Systems Biology of Drug Metabolizing Enzymes and Transporters: Relevance to Quantitative Systems Pharmacology.* Clin Pharmacol Ther 2020. [PMID 32119114](https://pubmed.ncbi.nlm.nih.gov/32119114/)
131. Elmokadem A, et al. *Quantitative Systems Pharmacology and Physiologically-Based Pharmacokinetic Modeling With mrgsolve: A Hands-On Tutorial.* CPT Pharmacometrics Syst Pharmacol 2019. [PMID 31652028](https://pubmed.ncbi.nlm.nih.gov/31652028/)
132. Lu T, et al. *gPKPDviz: A flexible R shiny tool for pharmacokinetic/pharmacodynamic simulations using mrgsolve.* CPT Pharmacometrics Syst Pharmacol 2024. [PMID 38082557](https://pubmed.ncbi.nlm.nih.gov/38082557/)
133. Ghallab A, et al. *Enteronephrohepatic circulation of bile acids and therapeutic potential of systemic bile acid transporter inhibition.* J Hepatol 2025. [PMID 40414504](https://pubmed.ncbi.nlm.nih.gov/40414504/)
134. De Bruijn VMP, et al. *Intestinal in vitro transport assay combined with physiologically based kinetic modeling.* ALTEX 2024. [PMID 37528756](https://pubmed.ncbi.nlm.nih.gov/37528756/)

---

## The model's claims mapped to their evidence

| Claim of the model | Type of evidence | References |
|---|---|---|
| Normal pool 3 g, 4–6 cycles per day, synthesis 0.35 g/day | 📖 calibration | 1, 7, 8, 10 |
| Ileal conservation f ≈ 0.97, ASBT distally weighted | 📖 calibration | 15–18, 102 |
| Colonic secretion threshold ≈ 3 mM dihydroxy bile acid | 📖 parameter taken directly | **52** (Mekjian 1971) |
| CA is weak and DCA the strongest secretagogue | 📖 potency weights | 52, 55, 60 |
| FGF19 is the main arm of CYP7A1 suppression, SHP auxiliary | 📖 structure | 27, 28, 29, 30 |
| Idiopathic (type 2) BAD = low fasting FGF19 (`κ↓`) | 📖 lesion axis | 66, 67, 70 |
| **At steady state the colonic load = hepatic synthesis** | ⚙ derived by the model | (mass conservation) |
| **Freeze the feedback and malabsorption alone gives no diarrhoea** | ⚙ derived by the model | to be tested |
| **SeHCAT reads `φ` only; C4/faecal bile acids read `S`** | ⚙ derived by the model | consistent with 69, 67, 68 |
| **SeHCAT-negative bile acid diarrhoea must exist** | ⚙ prediction | consistent with 68, 71 |
| A sequestrant raises C4 (structural escape) | ⚙ + 📖 | 77, 78 |
| **An FXR agonist raises SeHCAT without fixing the ileum** | ⚙ prediction (unverified) | 85, 86 |
| ASBT inhibitor = iatrogenic BAM | ⚙ + 📖 | 94–98 |
| **The 100 cm rule is derived from three constants** | ⚙ derived by the model | consistent with 102, 103, 104 |
| Loperamide reduces faecal bile acids too | ⚙ prediction | consistent with 62, 120, 121 |
| The microbiota changes potency, not load | ⚙ prediction | 109–113, 116 |

---

*134 papers · 134 PubMed links · every PMID looked up and verified through the NCBI E-utilities.*
