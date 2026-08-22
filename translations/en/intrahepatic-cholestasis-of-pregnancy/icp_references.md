# Intrahepatic cholestasis of pregnancy (ICP) — QSP model references

> **How this list was built.** Every PMID below was actually looked up through NCBI E-utilities (`esearch` +
> `esummary`), and **only the records that came back** have been transcribed here. First author, year
> and journal name are exactly as PubMed returned them; no entry was written from
> memory or from guesswork. The `→` after each entry states **which parameter or which
> structural decision** in the model that reference was used for. Values with no supporting evidence are
> stated explicitly as a "model choice".
>
> **How this list was built.** Every PMID below was resolved against NCBI
> E-utilities and only records that came back are listed; first author, year
> and journal are as PubMed returned them. The `→` note on each entry says
> which parameter or which structural decision in the model it supports.
> Where the model has no literature anchor, that is stated instead of
> concealed.

---

## 1. Diagnostic criteria · practice guidelines · epidemiology

The reference points for the model's severity stratification (≥10 diagnosis · ≥40 severe · ≥100 most severe) and for the
delivery-timing comparison.

1. Dajti E, 2025, *Nat Rev Dis Primers* — Intrahepatic cholestasis of pregnancy. [PMID 40707479](https://pubmed.ncbi.nlm.nih.gov/40707479/) → the most recent synthesis of the disease as a whole; the reference account of diagnostic thresholds, epidemiology and the current state of treatment
2. European Association for the Study of the Liver, 2023, *J Hepatol* — EASL Clinical Practice Guidelines on the management of liver diseases in pregnancy. [PMID 37394016](https://pubmed.ncbi.nlm.nih.gov/37394016/) → diagnostic criteria and the delivery-timing recommendation (the comparator for section G of the model)
3. Society for Maternal-Fetal Medicine (SMFM), 2021, *Am J Obstet Gynecol* — Consult Series #56: Intrahepatic cholestasis of pregnancy. [PMID 33197417](https://pubmed.ncbi.nlm.nih.gov/33197417/) → recommends delivery at 36–37 weeks above ≥100 µmol/L; contrasted with the model's optimisation result
4. Hobson SR, 2024, *J Obstet Gynaecol Can* — Guideline No. 452: Diagnosis and Management of Intrahepatic Cholestasis of Pregnancy. [PMID 39089469](https://pubmed.ncbi.nlm.nih.gov/39089469/) → the Canadian guideline, management by stratum
5. Hague WM, 2023, *Aust N Z J Obstet Gynaecol* — Diagnosis and management: consensus statement (SOMANZ). [PMID 37431680](https://pubmed.ncbi.nlm.nih.gov/37431680/) → stratification criteria
6. Kothari S, 2024, *Gastroenterology* — AGA Clinical Practice Update on Pregnancy-Related Gastrointestinal and Liver Disease. [PMID 39140906](https://pubmed.ncbi.nlm.nih.gov/39140906/) → management recommendations from a gastroenterological perspective
7. Giouleka S, 2025, *Obstet Gynecol Surv* — Intrahepatic Cholestasis of Pregnancy: A Comparative Review of Guidelines. [PMID 41194372](https://pubmed.ncbi.nlm.nih.gov/41194372/) → the extent of disagreement between guidelines (the grounds for not fitting the model to any single guideline)
8. Bicocca MJ, 2018, *Eur J Obstet Gynecol Reprod Biol* — Review of six national and regional guidelines. [PMID 30396107](https://pubmed.ncbi.nlm.nih.gov/30396107/) → same point
9. Mitchell AL, 2021, *BJOG* — Re-evaluating diagnostic thresholds for intrahepatic cholestasis of pregnancy. [PMID 33586324](https://pubmed.ncbi.nlm.nih.gov/33586324/) → a non-fasting threshold of ≥19 µmol/L; why the model reports the "assay value" and the "endogenous fraction" separately
10. Rodriguez N, 2026, *Am J Gastroenterol* — Evaluating the 2023 EASL Guidelines: Risk Stratification and Outcomes. [PMID 41351236](https://pubmed.ncbi.nlm.nih.gov/41351236/) → external validation of the stratification
11. Tran TT, 2016, *Am J Gastroenterol* — ACG Clinical Guideline: Liver Disease and Pregnancy. [PMID 26832651](https://pubmed.ncbi.nlm.nih.gov/26832651/) → placement within the wider context of liver disease
12. Contreras Vidal C, 2026, *Medwave* — ICP in Chile: epidemiological change and a microbiological hypothesis. [PMID 42008798](https://pubmed.ncbi.nlm.nih.gov/42008798/) → variation in incidence between regions and eras (represented in the model by the genetic susceptibility vector)
13. Piechota J, 2020, *J Clin Med* — Intrahepatic Cholestasis in Pregnancy: Review of the Literature. [PMID 32384779](https://pubmed.ncbi.nlm.nih.gov/32384779/) → background account

---

## 2. Perinatal outcomes and stillbirth risk — the three numbers the model was fitted to

The first reference in this section is the source of the **only three constants in the model fitted to outcome epidemiology**
(`HSB0`, `HSBSC`, `HN`).

14. **Ovadia C, 2019, *Lancet*** — Association of adverse perinatal outcomes of ICP with biochemical markers: aggregate and individual patient data meta-analyses. [PMID 30773280](https://pubmed.ncbi.nlm.nih.gov/30773280/) → **stillbirth 0.13% / 0.28% / 3.44% (TBA <40 / 40–99 / ≥100 µmol/L).** Section C of `icp_calibration.py` fits `HSB0`·`HSBSC`·`HN` to these three values, and section D refits the same three values against three different driving variables to measure where the steep part of the threshold lies
15. Sarker M, 2022, *Am J Obstet Gynecol* — Beyond stillbirth: association of ICP severity and adverse outcomes. [PMID 36008054](https://pubmed.ncbi.nlm.nih.gov/36008054/) → frequencies by stratum of the endpoints other than stillbirth (meconium · preterm birth · NICU) → the reference for `HM0`, `HMSC`, `HPT0`
16. Zhou Q, 2024, *BMC Pregnancy Childbirth* — Severity of ICP increases risks of adverse outcomes beyond stillbirth: 15,826 patients. [PMID 38997626](https://pubmed.ncbi.nlm.nih.gov/38997626/) → same purpose, large-scale validation
17. Huang X, 2024, *PLoS One* — Systematic review and meta-analysis of obstetric and neonatal outcomes. [PMID 38833446](https://pubmed.ncbi.nlm.nih.gov/38833446/) → the reference range for the NICU admission rate
18. Yao L, 2025, *BMC Pregnancy Childbirth* — Total bile acid concentrations and adverse perinatal outcomes, including asymptomatic hypercholanemia. [PMID 41102717](https://pubmed.ncbi.nlm.nih.gov/41102717/) → the existence of asymptomatic hypercholanaemia; the clinical grounds for the itch axis and the bile acid axis having to be separate in the model
19. Çalık MG, 2026, *BMC Pregnancy Childbirth* — Perinatal outcomes in relation to transaminase and bile acid levels. [PMID 41987074](https://pubmed.ncbi.nlm.nih.gov/41987074/) → the observation that ALT and TBA carry different information
20. Granese R, 2023, *J Clin Med* — Maternal and Neonatal Outcomes in ICP. [PMID 37445442](https://pubmed.ncbi.nlm.nih.gov/37445442/) → the frequency of meconium staining
21. Feng F, 2024, *Sci Rep* — Clinical subtypes and bile acid levels with pregnancy outcomes. [PMID 38806569](https://pubmed.ncbi.nlm.nih.gov/38806569/) → outcomes by subtype
22. Sarker MR, 2025, *Am J Perinatol* — Adverse Outcomes Associated with Progressive ICP. [PMID 39592109](https://pubmed.ncbi.nlm.nih.gov/39592109/) → outcomes where the disease worsens over time; directly relevant to the "smooth trajectory" limitation noted in section G of the model
23. Sarker MR, 2025, *Pregnancy (Hoboken)* — Earlier diagnosis of ICP and adverse pregnancy outcomes. [PMID 40959761](https://pubmed.ncbi.nlm.nih.gov/40959761/) → the prognosis of early presentation
24. Hocaoglu M, 2025, *Fetal Pediatr Pathol* — Early- and Late-Onset ICP: maternal and neonatal outcomes. [PMID 39932778](https://pubmed.ncbi.nlm.nih.gov/39932778/) → time of onset and outcome
25. Niculae LE, 2025, *Biomedicines* — Neonatal Impact Through the Lens of Current Evidence. [PMID 41007629](https://pubmed.ncbi.nlm.nih.gov/41007629/) → neonatal outcomes
26. Nielsen JH, 2021, *Acta Obstet Gynecol Scand* — Differentiated timing of induction for women with ICP. [PMID 32970824](https://pubmed.ncbi.nlm.nih.gov/32970824/) → observational data on the timing of induction by stratum; the comparator for section G of the model
27. Zehner L, 2023, *Arch Gynecol Obstet* — Evaluation of obstetric management in German maternity units. [PMID 36030428](https://pubmed.ncbi.nlm.nih.gov/36030428/) → the distribution of delivery timing in actual practice
28. Xu T, 2024, *Int J Gynaecol Obstet* — Risk-stratified management strategies: tertiary centre review. [PMID 37470272](https://pubmed.ncbi.nlm.nih.gov/37470272/) → stratified management in practice
29. Capatina N, 2024, *Obstet Med* — Meta-analyses in cholestatic pregnancy: the outstanding clinical questions. [PMID 39262915](https://pubmed.ncbi.nlm.nih.gov/39262915/) → the list of unresolved questions; the statement of the two paradoxes this model targets

---

## 3. Genetic susceptibility

The evidence behind the model's five patient descriptors `GBSEP`, `GMDR3`, `GFIC1`, `GSULT`, `GFXR`.

30. Dixon PH, 2022, *Nat Commun* — GWAS meta-analysis of ICP implicates multiple hepatic genes and regulatory elements. [PMID 35977952](https://pubmed.ncbi.nlm.nih.gov/35977952/) → that ICP is polygenic; why the model uses a continuous multiplier on transport capacity rather than a single variant
31. Tyrmi JS, 2026, *Nat Commun* — Genome-wide meta-analysis identifies genetic drivers of bile acid metabolism in ICP. [PMID 42178310](https://pubmed.ncbi.nlm.nih.gov/42178310/) → same point, centred on the bile acid metabolism axis
32. Dixon PH, 2009, *Gut* — Contribution of variant alleles of ABCB11 to susceptibility to ICP. [PMID 18987030](https://pubmed.ncbi.nlm.nih.gov/18987030/) → the evidence for **`GBSEP` (reduced ABCB11/BSEP transport capacity)**
33. Zhang D, 2025, *Sci Rep* — Association between ABCB4 variants and intrahepatic cholestasis of pregnancy. [PMID 39865141](https://pubmed.ncbi.nlm.nih.gov/39865141/) → the evidence for **`GMDR3`**
34. Prescher M, 2019, *Biol Chem* — ABCB4/MDR3 in health and disease. [PMID 30730833](https://pubmed.ncbi.nlm.nih.gov/30730833/) → the functional evidence that MDR3 is a phospholipid floppase → the structure of the model's biliary bile acid:PC ratio (`BSPC`)
35. Sticova E, 2020, *Ann Hepatol* — ABCB4 disease: many faces of one gene deficiency. [PMID 31759867](https://pubmed.ncbi.nlm.nih.gov/31759867/) → the ABCB4 phenotypic spectrum (including whether GGT is raised)
36. Nayagam JS, 2020, *Aliment Pharmacol Ther* — Liver disease in adults with variants in ABCB11, ABCB4 and ATP8B1. [PMID 33070363](https://pubmed.ncbi.nlm.nih.gov/33070363/) → the adult phenotypes of the three genes, **`GFIC1` (ATP8B1)** included
37. Yeap SP, 2019, *J Gastroenterol Hepatol* — Biliary transporter gene mutations in severe ICP. [PMID 29992621](https://pubmed.ncbi.nlm.nih.gov/29992621/) → the observation that variant frequency is higher in severe ICP; the model's severity-genotype correspondence
38. Aydın GA, 2020, *Taiwan J Obstet Gynecol* — The role of genetic mutations in ICP. [PMID 32917322](https://pubmed.ncbi.nlm.nih.gov/32917322/) → same point
39. Stieger B, 2011, *Expert Opin Drug Metab Toxicol* — Genetic variations of bile salt transporters as predisposing factors. [PMID 21320040](https://pubmed.ncbi.nlm.nih.gov/21320040/) → the susceptibility shared by drug-induced cholestasis and ICP
40. Pauli-Magnus C, 2010, *Semin Liver Dis* — Genetic determinants of drug-induced cholestasis and ICP. [PMID 20422497](https://pubmed.ncbi.nlm.nih.gov/20422497/) → same point
41. Degiorgio D, 2016, *J Gastroenterol* — ABCB4 mutations in adult patients with cholestatic liver disease. [PMID 26324191](https://pubmed.ncbi.nlm.nih.gov/26324191/) → residual function by variant; the basis for the range of `GMDR3` values
42. Falcão D, 2022, *Dig Liver Dis* — The wide phenotypic and genetic spectrum of ABCB4 deficiency. [PMID 34376370](https://pubmed.ncbi.nlm.nih.gov/34376370/) → same point
43. Henkel SA, 2019, *World J Hepatol* — Expanding etiology of progressive familial intrahepatic cholestasis. [PMID 31183005](https://pubmed.ncbi.nlm.nih.gov/31183005/) → the PFIC spectrum (TJP2, NR1H4 included)

---

## 4. Sex hormones inhibit BSEP — the lesion itself

The two inhibition terms that generate the disease in the model (`KI_E2G`, `KI_P4S`), and why third-trimester
onset and resolution after delivery emerge **with no change of parameters**.

44. **Vallejo M, 2006, *J Hepatol*** — Potential role of trans-inhibition of the bile salt export pump by progesterone metabolites in the etiopathogenesis of ICP. [PMID 16458994](https://pubmed.ncbi.nlm.nih.gov/16458994/) → the direct evidence for **`KI_P4S` (trans-inhibition of BSEP by progesterone sulfates)**
45. **Abu-Hayyeh S, 2013, *Hepatology*** — ICP levels of sulfated progesterone metabolites inhibit FXR resulting in a cholestatic phenotype. [PMID 22961653](https://pubmed.ncbi.nlm.nih.gov/22961653/) → FXR antagonism by P4-S (the `P4S → FXR` inhibition edge on the model map)
46. Abu-Hayyeh S, 2010, *J Biol Chem* — Inhibition of NTCP-mediated bile acid transport by cholestatic sulfated progesterone metabolites. [PMID 20177056](https://pubmed.ncbi.nlm.nih.gov/20177056/) → NTCP inhibition; consistent with the impaired uptake represented by the model's `updep` term
47. Abu-Hayyeh S, 2016, *Hepatology* — Prognostic and mechanistic potential of progesterone sulfates in ICP and pruritus gravidarum. [PMID 26426865](https://pubmed.ncbi.nlm.nih.gov/26426865/) → the observation that **P4-S tracks itch better than bile acids do** → the key evidence behind the structural decision to drive autotaxin from P4S in the model (`EA_P4S`, `KA_P4S`, `HP4S`)
48. Meng LJ, 1997, *Hepatology* — Progesterone metabolites and bile acids in serum of patients with ICP: effect of UDCA therapy. [PMID 9398000](https://pubmed.ncbi.nlm.nih.gov/9398000/) → the observation that **UDCA fails to normalise progesterone sulfates**; why UDCA barely reduces itch in the model
49. Meng LJ, 1997, *J Hepatol* — Effects of UDCA on conjugated bile acids and progesterone metabolites in serum and urine of ICP. [PMID 9453429](https://pubmed.ncbi.nlm.nih.gov/9453429/) → same point
50. Meng LJ, 1997, *J Hepatol* — Profiles of bile acids and progesterone metabolites in urine and serum of women with ICP. [PMID 9288610](https://pubmed.ncbi.nlm.nih.gov/9288610/) → the P4-S concentration range → the magnitude of `KFP4S`, `KOUTP4S`
51. Reyes H, 2000, *Ann Med* — Bile acids and progesterone metabolites in ICP. [PMID 10766400](https://pubmed.ncbi.nlm.nih.gov/10766400/) → the parallel rise of the two axes
52. Mitchell AL, 2025, *Am J Physiol Gastrointest Liver Physiol* — Progesterone sulfates are enterohepatically recycled and stimulate GPBAR1-mediated gut hormone release. [PMID 39888313](https://pubmed.ncbi.nlm.nih.gov/39888313/) → that P4-S is enterohepatically recycled (simplified in the model to a plasma state variable only — **an explicit simplification**)
53. Sanchon-Sanchez P, 2024, *Biochim Biophys Acta Mol Basis Dis* — Relationship between cholestasis and altered progesterone metabolism in the placenta-maternal liver tandem. [PMID 37956602](https://pubmed.ncbi.nlm.nih.gov/37956602/) → progesterone metabolism in the placenta-maternal liver axis
54. Ren L, 2019, *Biochem Biophys Res Commun* — Epiallopregnanolone sulfate induces hepatic bile acid accumulation and liver injury. [PMID 31575408](https://pubmed.ncbi.nlm.nih.gov/31575408/) → experimental evidence that P4-S is causal
55. Liu H, 2023, *Nat Commun* — Structural basis of bile salt extrusion and small-molecule inhibition in human BSEP. [PMID 37949847](https://pubmed.ncbi.nlm.nih.gov/37949847/) → the structural basis of BSEP inhibition; the plausibility of the competitive-inhibition form
56. Hosamani S, 2024, *J Phys Chem Lett* — Cholesterol allosterically modulates the structure and dynamics of ABCB11. [PMID 39058973](https://pubmed.ncbi.nlm.nih.gov/39058973/) → evidence that the membrane environment alters BSEP activity → the structure in which `GFIC1` (ATP8B1) multiplies BSEP

---

## 5. Bile acid synthesis · feedback · transport

57. Chiang JY, 2009, *J Lipid Res* — Bile acids: regulation of synthesis. [PMID 19346330](https://pubmed.ncbi.nlm.nih.gov/19346330/) → `VSYN`, `FRCA` (the CYP7A1/CYP8B1 branch point) and the whole FXR-SHP feedback structure
58. Holt JA, 2003, *Genes Dev* — Definition of a novel growth factor-dependent signal cascade for the suppression of bile acid biosynthesis. [PMID 12815072](https://pubmed.ncbi.nlm.nih.gov/12815072/) → the ileal FXR–FGF19 axis → `KS_F19`, `EC50FXI`, `KF19`
59. Song KH, 2009, *Hepatology* — Bile acids activate FGF19 signaling in human hepatocytes to inhibit CYP7A1. [PMID 19085950](https://pubmed.ncbi.nlm.nih.gov/19085950/) → FGF19–CYP7A1 inhibition in human cells
60. Zollner G, 2007, *Liver Int* — Expression of bile acid synthesis and detoxification enzymes and the alternative efflux pump MRP4 in PBC. [PMID 17696930](https://pubmed.ncbi.nlm.nih.gov/17696930/) → human data showing that **MRP4 is induced in cholestasis** → `EM4`, `KM4`, `HM4` (the evidence for a Hill form that is quiet at baseline and strongly induced in cholestasis)
61. Schaffner CA, 2015, *Liver Int* — The organic solute transporters alpha and beta are induced by hypoxia in human hepatocytes. [PMID 24703425](https://pubmed.ncbi.nlm.nih.gov/24703425/) → the inducibility of OSTα-β
62. Xu S, 2014, *Am J Physiol Gastrointest Liver Physiol* — A novel RARα/CAR-mediated mechanism for regulation of human OSTβ. [PMID 24264050](https://pubmed.ncbi.nlm.nih.gov/24264050/) → transcriptional regulation of the basolateral efflux route
63. van Dijk R, 2015, *Liver Int* — Characterization and treatment of persistent hepatocellular secretory failure. [PMID 24905729](https://pubmed.ncbi.nlm.nih.gov/24905729/) → the clinical phenotype of hepatocellular secretory failure
64. Kastrinou Lampou V, 2023, *Toxicol In Vitro* — Novel insights into bile acid detoxification via CYP, UGT and SULT enzymes. [PMID 36473578](https://pubmed.ncbi.nlm.nih.gov/36473578/) → the evidence for the two detoxification routes **`V6A` (CYP3A4 6α-hydroxylation)** and **`VLCAS` (SULT2A1 LCA sulfation)**
65. Bansal S, 2016, *J Chromatogr B* — Fast and sensitive quantification of human liver cytosolic lithocholic acid sulfation. [PMID 26773894](https://pubmed.ncbi.nlm.nih.gov/26773894/) → human hepatic LCA sulfation capacity → the magnitude of `VLCAS`, `KMLCAS`
66. Mao F, 2019, *J Biol Chem* — Increased sulfation of bile acids in NTCP deficiency. [PMID 31201272](https://pubmed.ncbi.nlm.nih.gov/31201272/) → when bile acid excretion is blocked, sulfation and renal excretion become the alternative route → `ESULT`
67. Chatterjee B, 2005, *Methods Enzymol* — Vitamin D receptor regulation of SULT2A1. [PMID 16399349](https://pubmed.ncbi.nlm.nih.gov/16399349/) → LCA–VDR–SULT2A1 feedback (shown on the map, folded into the PXR term in the ODEs — **an explicit simplification**)
68. Mok HY, 1977, *Gastroenterology* — Regulation of pool size of bile acids in man. [PMID 892372](https://pubmed.ncbi.nlm.nih.gov/892372/) → **the human bile acid pool size and the number of cycles per day** → sets the magnitude of `KBILE`, `KTR`, `VASBT`
69. Low-Beer TS, 1973, *Br Med J* — Regulation of bile salt pool size in man. [PMID 4704519](https://pubmed.ncbi.nlm.nih.gov/4704519/) → same purpose
70. Hepner GW, 1975, *Gastroenterology* — Effect of decreased gallbladder stimulation on enterohepatic cycling and kinetics of bile acids. [PMID 1132637](https://pubmed.ncbi.nlm.nih.gov/1132637/) → gallbladder storage and cycling frequency
71. Angelin B, 1978, *J Lipid Res* — Individual serum bile acid concentrations in normo- and hyperlipoproteinemia: relation to pool size. [PMID 670831](https://pubmed.ncbi.nlm.nih.gov/670831/) → **the normal plasma bile acid concentration (the reference for the model's 4–5 µmol/L in normal pregnancy)**
72. Hofmann AF, 1987, *Gastroenterology* — Simulation of the metabolism and enterohepatic circulation of endogenous deoxycholic acid in humans using a physiologic pharmacokinetic model. [PMID 3623017](https://pubmed.ncbi.nlm.nih.gov/3623017/) → the direct antecedent of this model's enterohepatic compartment structure; the treatment of colonic 7α-dehydroxylation (`KDH_CA`)
73. Cravetto C, 1988, *Hepatology* — Computer simulation of portal venous shunting and other isolated hepatobiliary defects of the enterohepatic circulation using a physiologic model. [PMID 3391514](https://pubmed.ncbi.nlm.nih.gov/3391514/) → the modelling precedent for defects of hepatic first-pass extraction → `EUP*`, `JHALF`, `KHLSAT`
74. Hulzebos CV, 2001, *J Lipid Res* — Measurement of parameters of cholic acid kinetics in plasma using stable isotope dilution. [PMID 11714862](https://pubmed.ncbi.nlm.nih.gov/11714862/) → cholic acid kinetic parameters
75. Rutgeerts P, 1983, *J Lipid Res* — The enterohepatic circulation of bile acids during continuous liquid formula perfusion of the duodenum. [PMID 6875385](https://pubmed.ncbi.nlm.nih.gov/6875385/) → intraluminal bile acid concentration (of the order of mmol/L) → `MIC50` and the micelle formation term

---

## 6. Placental bile acid transport

The model's `VP`, `KMP*`, `PS*`, `TRANSIN`. **The references in this section support placental transport being
saturable and bidirectional, but no study directly measuring the absolute value of `VP` in humans
was found — `VP` is a model choice fitted to the cord blood concentration of normal
pregnancy.**

76. Monte MJ, 1995, *Pediatr Res* — Relationship between bile acid transplacental gradients and transport across the fetal-facing plasma membrane of the human trophoblast. [PMID 7478809](https://pubmed.ncbi.nlm.nih.gov/7478809/) → **the relationship between the fetal-maternal concentration gradient and trophoblast transport** → the structure of `PS*` and `JF2M`
77. Marin JJ, 1995, *Am J Physiol* — ATP-dependent bile acid transport across microvillous membrane of human term trophoblast. [PMID 7733292](https://pubmed.ncbi.nlm.nih.gov/7733292/) → the existence of active (ATP-dependent) transport → a saturable `VP`
78. Marín JJ, 2005, *Ann Hepatol* — Molecular bases of the excretion of fetal bile acids and pigments through the fetal liver-placenta-maternal liver pathway. [PMID 16010240](https://pubmed.ncbi.nlm.nih.gov/16010240/) → **the structure in which the placenta is the only substantive excretion route for fetal bile acids** → the design of the escape route out of the fetal compartment in the model
79. Macias RI, 2000, *Hepatology* — Effect of maternal cholestasis on bile acid transfer across the rat placenta-maternal liver tandem. [PMID 10733555](https://pubmed.ncbi.nlm.nih.gov/10733555/) → maternal cholestasis impairs placental transfer → the `CAPUTIL` saturation concept
80. Macias RI, 2009, *World J Gastroenterol* — Excretion of biliary compounds during intrauterine life. [PMID 19230042](https://pubmed.ncbi.nlm.nih.gov/19230042/) → a synthesis of the fetal excretion routes
81. Ontsouka E, 2021, *Int J Mol Sci* — Placental Expression of Bile Acid Transporters in Intrahepatic Cholestasis of Pregnancy. [PMID 34638773](https://pubmed.ncbi.nlm.nih.gov/34638773/) → altered transporter expression in the ICP placenta (human) → the grounds for making `GP` a patient descriptor
82. St-Pierre MV, 2000, *Am J Physiol Regul Integr Comp Physiol* — Expression of members of the multidrug resistance protein family in human term placenta. [PMID 11004020](https://pubmed.ncbi.nlm.nih.gov/11004020/) → placental MRP expression
83. Blazquez AG, 2014, *Toxicol Appl Pharmacol* — Acetaminophen impairs the placental barrier to bile acids via BCRP. [PMID 24631341](https://pubmed.ncbi.nlm.nih.gov/24631341/) → evidence that the placental barrier can be modulated by drugs (the model's `GP` manipulation scenario)
84. Serrano MA, 1998, *J Hepatol* — Beneficial effect of UDCA on alterations induced by cholestasis of pregnancy in bile acid transport across the human placenta. [PMID 9625319](https://pubmed.ncbi.nlm.nih.gov/9625319/) → **direct evidence that UDCA intervenes in placental transport** → the source of `KMP5`, `TRANSIN` (competition by UDCA for the transporter)
85. Serrano MA, 2003, *J Pharmacol Exp Ther* — Effect of UDCA on the impairment induced by maternal cholestasis in the rat placenta-maternal liver tandem. [PMID 12606635](https://pubmed.ncbi.nlm.nih.gov/12606635/) → same point
86. Estiú MC, 2015, *Br J Clin Pharmacol* — Effect of UDCA on the altered progesterone and bile acid homeostasis in the mother-placenta-foetus trio. [PMID 25099365](https://pubmed.ncbi.nlm.nih.gov/25099365/) → the effects of UDCA by bile acid species in the mother-placenta-fetus trio
87. Basu S, 2025, *Am J Physiol Gastrointest Liver Physiol* — Unresolved alterations in bile acid composition and dyslipidemia in maternal and cord blood after UDCA treatment for ICP. [PMID 39947696](https://pubmed.ncbi.nlm.nih.gov/39947696/) → **the abnormal cord blood bile acid composition persists even after UDCA treatment** → direct experimental support for the part of the model that says "the assay value comes down but the fetal hydrophobic burden does not come down as much"
88. Wang G, 2026, *Mol Cell Biochem* — Bile acid-induced TBG depletion promotes trophoblast apoptosis in ICP. [PMID 42371391](https://pubmed.ncbi.nlm.nih.gov/42371391/) → trophoblast apoptosis (`TROPH`)
89. Wu WB, 2015, *Placenta* — FXR agonist protects against bile acid induced damage and oxidative stress in mouse placenta. [PMID 25747729](https://pubmed.ncbi.nlm.nih.gov/25747729/) → placental oxidative stress
90. Keitel V, 2013, *Placenta* — Effect of maternal cholestasis on TGR5 expression in human and rat placenta at term. [PMID 23849932](https://pubmed.ncbi.nlm.nih.gov/23849932/) → placental TGR5

---

## 7. The fetal myocardium — where the threshold actually is

**The most important section of the model.** `IMAXGJ`, `IC50GJ`, `HGJ`, `ECA`, `KCA` all
come from this section, and the ablation experiment in section I of `icp_calibration.py` shows that **replacing the Hill form
here with a linear one makes the 100 µmol/L threshold disappear**.

91. **Williamson C, 2001, *Clin Sci (Lond)*** — The bile acid taurocholate impairs rat cardiomyocyte function: a proposed mechanism for intra-uterine fetal death in obstetric cholestasis. [PMID 11256973](https://pubmed.ncbi.nlm.nih.gov/11256973/) → **the starting point of the entire fetal myocardium module of the model**
92. **Gorelik J, 2002, *Clin Sci (Lond)*** — Taurocholate induces changes in rat cardiomyocyte contraction and calcium dynamics. [PMID 12149111](https://pubmed.ncbi.nlm.nih.gov/12149111/) → the direct evidence for **`ECA`, `KCA` (calcium overload)**
93. **Gorelik J, 2004, *BJOG*** — Comparison of the arrhythmogenic effects of tauro- and glycoconjugates of cholic acid in an in vitro study of rat cardiomyocytes. [PMID 15270939](https://pubmed.ncbi.nlm.nih.gov/15270939/) → **arrhythmogenic potency differs between bile acid species** → the grounds for the model's species-specific cytotoxicity weights `WTOX`
94. Sheikh Abdul Kadir SH, 2010, *PLoS One* — Bile acid-induced arrhythmia is mediated by muscarinic M2 receptors in neonatal rat cardiomyocytes. [PMID 20300620](https://pubmed.ncbi.nlm.nih.gov/20300620/) → the M2 receptor route (`MRECEPT` on the map)
95. Ibrahim E, 2018, *Sci Rep* — Bile acids and their respective conjugates elicit different responses in neonatal cardiomyocytes: role of Gi protein, muscarinic receptors and TGR5. [PMID 29740092](https://pubmed.ncbi.nlm.nih.gov/29740092/) → quantitative evidence for differences in response by species and by conjugate → the ordering of `WTOX`
96. Miragoli M, 2011, *Hepatology* — A protective antiarrhythmic role of ursodeoxycholic acid in an in vitro rat model of the cholestatic fetal heart. [PMID 21809354](https://pubmed.ncbi.nlm.nih.gov/21809354/) → **the direct antiarrhythmic effect of UDCA**; the model represents this only as `WTOX[UDCA] = 0.02` and includes no separate myocardial protection term — **stated explicitly as a point at which the model may underestimate the fetal cardiac protection given by UDCA**
97. Adeyemi O, 2017, *PLoS One* — UDCA prevents ventricular conduction slowing and arrhythmia by restoring T-type calcium current in fetuses during cholestasis. [PMID 28934223](https://pubmed.ncbi.nlm.nih.gov/28934223/) → same point, in vivo fetus
98. Gorelik J, 2003, *BJOG* — Dexamethasone and UDCA protect against the arrhythmogenic effect of taurocholate. [PMID 12742331](https://pubmed.ncbi.nlm.nih.gov/12742331/) → in vitro protection by dexamethasone (the clinical RCT was negative — marked as such on the map)
99. Gorelik J, 2006, *BJOG* — Genes encoding bile acid, phospholipid and anion transporters are expressed in a human fetal cardiomyocyte culture. [PMID 16637898](https://pubmed.ncbi.nlm.nih.gov/16637898/) → fetal cardiomyocytes express bile acid transporters (human)
100. Abdul Kadir SH, 2009, *J Cell Mol Med* — Embryonic stem cell-derived cardiomyocytes as a model to study fetal arrhythmia related to maternal disease. [PMID 19438812](https://pubmed.ncbi.nlm.nih.gov/19438812/) → a human-derived model
101. Schultz F, 2016, *Prog Biophys Mol Biol* — The protective effect of UDCA in an in vitro model of the human fetal heart occurs via targeting cardiac fibroblasts. [PMID 26777584](https://pubmed.ncbi.nlm.nih.gov/26777584/) → a fibroblast-mediated mechanism (not included in the model — **an explicit omission**)
102. Williamson C, 2011, *Dig Dis* — Bile acid signaling in fetal tissues: implications for ICP. [PMID 21691106](https://pubmed.ncbi.nlm.nih.gov/21691106/) → a synthesis of bile acid signalling in fetal tissues
103. Guerra M, 2025, *Am J Perinatol* — Fetal TEI Index in Pregnancies with ICP: A Case-Control Study. [PMID 40112872](https://pubmed.ncbi.nlm.nih.gov/40112872/) → **the clinical observation that indices of cardiac function really do change in the human fetus** → the grounds for holding that the model's myocardial module is testable in humans as well
104. Şık ME, 2026, *J Perinat Med* — Evaluation of Doppler parameters and obstetric outcomes in ICP. [PMID 41770037](https://pubmed.ncbi.nlm.nih.gov/41770037/) → the observation that Doppler does not predict acute events (consistent with the clinical implication of the model)
105. Sepúlveda WH, 1991, *Eur J Obstet Gynecol Reprod Biol* — Vasoconstrictive effect of bile acids on isolated human placental chorionic veins. [PMID 1773876](https://pubmed.ncbi.nlm.nih.gov/1773876/) → the direct evidence for **`EV`, `KV` (chorionic plate vein constriction)**

---

## 8. The itch axis — it is separate from the bile acid axis

The model's second axis (`ATX`, `LPA`, `ITCHC`). **The central observation of this section is that "plasma
autotaxin correlates with itch but total bile acids do not", and that observation is
why `EA_CH` is set small and `KA_CH` low (that is, saturated) in the model.**

106. **Kremer AE, 2010, *Gastroenterology*** — Lysophosphatidic acid is a potential mediator of cholestatic pruritus. [PMID 20546739](https://pubmed.ncbi.nlm.nih.gov/20546739/) → the first proposal that LPA is the mediator → the model's `LPA` state variable
107. **Kremer AE, 2012, *Hepatology*** — Serum autotaxin is increased in pruritus of cholestasis, but not of other origin, and responds to therapeutic interventions. [PMID 22473838](https://pubmed.ncbi.nlm.nih.gov/22473838/) → **the source of `ERIFA` (rifampicin lowers autotaxin) and of the observation that UDCA does not.** The entire two-axis separation structure of the model hangs on this paper
108. Beuers U, 2023, *Nat Rev Gastroenterol Hepatol* — Mechanisms of pruritus in cholestasis: understanding and treating the itch. [PMID 36307649](https://pubmed.ncbi.nlm.nih.gov/36307649/) → the most recent synthesis; the axis on which each treatment acts
109. Beuers U, 2014, *Hepatology* — Pruritus in cholestasis: facts and fiction. [PMID 24807046](https://pubmed.ncbi.nlm.nih.gov/24807046/) → states explicitly how weak the bile acid-itch correlation is
110. Oude Elferink RP, 2011, *Dig Dis* — The molecular mechanism of cholestatic pruritus. [PMID 21691108](https://pubmed.ncbi.nlm.nih.gov/21691108/) → a synthesis of the mechanism
111. Oude Elferink RP, 2015, *Biochim Biophys Acta* — Lysophosphatidic acid and signaling in sensory neurons. [PMID 25218302](https://pubmed.ncbi.nlm.nih.gov/25218302/) → LPA→sensory neurone → the evidence for `HLP` (the apparent cooperativity of a multi-step amplification)
112. Wunsch E, 2016, *Sci Rep* — Serum autotaxin is a marker of the severity of liver injury and overall survival in cholestatic liver diseases. [PMID 27506882](https://pubmed.ncbi.nlm.nih.gov/27506882/) → the autotaxin concentration range (the model's ATX multiplier of 1.5–4)
113. Keune WJ, 2016, *Nat Commun* — Steroid binding to autotaxin links bile salts and lysophosphatidic acid signalling. [PMID 27075612](https://pubmed.ncbi.nlm.nih.gov/27075612/) → **bile acids bind autotaxin directly** → the structural basis for the model's weak `EA_CH` term
114. Langedijk JAGM, 2021, *Biochim Biophys Acta Mol Basis Dis* — Inhibition of autotaxin by bile salts and bile salt-like molecules increases its expression by feedback regulation. [PMID 34389475](https://pubmed.ncbi.nlm.nih.gov/34389475/) → **when bile acids inhibit autotaxin, expression rises through feedback** → why the bile acid-itch relationship is not monotonic; consistent with the model's saturating term
115. **Yu H, 2019, *Elife*** — MRGPRX4 is a bile acid receptor for human cholestatic itch. [PMID 31500698](https://pubmed.ncbi.nlm.nih.gov/31500698/) → the evidence for **`WBA`, `KBA` (the direct bile acid-itch route)**
116. Meixiong J, 2019, *Proc Natl Acad Sci U S A* — MRGPRX4 is a GPCR activated by bile acids that may contribute to cholestatic pruritus. [PMID 31068464](https://pubmed.ncbi.nlm.nih.gov/31068464/) → an independent report of the same finding
117. Yu H, 2021, *Semin Liver Dis* — MRGPRX4 in Cholestatic Pruritus. [PMID 34161994](https://pubmed.ncbi.nlm.nih.gov/34161994/) → a synthesis
118. Rodrigo M, 2026, *Liver Int* — Defining the contribution of genetic variants in MRGPRX4 with pruritus in paediatric cholestasis. [PMID 41817014](https://pubmed.ncbi.nlm.nih.gov/41817014/) → MRGPRX4 variants change the degree of itch (not included in the model — **stated as a candidate for extension**)
119. Yang J, 2026, *Nat Chem Biol* — Development of a clinically viable MRGPRX4 inverse agonist for cholestatic itch treatment. [PMID 41957282](https://pubmed.ncbi.nlm.nih.gov/41957282/) → a new drug aimed directly at the itch axis; the route the model's `WBA` term takes as its target
120. Yang J, 2024, *Cell* — Structure-guided discovery of bile acid derivatives for treating liver diseases without causing itch. [PMID 39476841](https://pubmed.ncbi.nlm.nih.gov/39476841/) → direct evidence that the bile acid axis and the itch axis are pharmacologically separable
121. Wolters F, 2025, *Biochim Biophys Acta Mol Basis Dis* — The role of bile salts in itch receptor activation. [PMID 40628364](https://pubmed.ncbi.nlm.nih.gov/40628364/) → receptor activation by bile acid species
122. Alemi F, 2013, *J Clin Invest* — The TGR5 receptor mediates bile acid-induced itch and analgesia. [PMID 23524965](https://pubmed.ncbi.nlm.nih.gov/23524965/) → the TGR5 route (shown on the map)
123. Lieu T, 2014, *Gastroenterology* — The bile acid receptor TGR5 activates the TRPA1 channel to induce itch in mice. [PMID 25194674](https://pubmed.ncbi.nlm.nih.gov/25194674/) → TRPA1 downstream
124. Song MH, 2022, *Biomol Ther (Seoul)* — Lithocholic acid activates MAS-related GPCRs, contributing to itch in mice. [PMID 34263729](https://pubmed.ncbi.nlm.nih.gov/34263729/) → the evidence that **LCA is the most strongly pruritogenic species**; the consistency between `WTOX` and `WBA`
125. Akiyama T, 2013, *Neuroscience* — Neural processing of itch. [PMID 23891755](https://pubmed.ncbi.nlm.nih.gov/23891755/) → spinal GRPR and central sensitisation → `KSENS`, `KSC` (the scratch-itch feedback)
126. Vander Does A, 2022, *Am J Clin Dermatol* — Cholestatic Itch: Our Current Understanding of Pathophysiology and Treatments. [PMID 35900649](https://pubmed.ncbi.nlm.nih.gov/35900649/) → the grounds for assigning each treatment to an axis
127. Kanda T, 2025, *Int J Mol Sci* — Pruritus in Chronic Cholestatic Liver Diseases. [PMID 40076514](https://pubmed.ncbi.nlm.nih.gov/40076514/) → the most recent synthesis
128. Trivella J, 2021, *Expert Opin Drug Saf* — Safety considerations for the management of cholestatic itch. [PMID 33836644](https://pubmed.ncbi.nlm.nih.gov/33836644/) → safety (the constraints in pregnancy)
129. Düll MM, 2019, *Curr Gastroenterol Rep* — Treatment of Pruritus Secondary to Liver Disease. [PMID 31367993](https://pubmed.ncbi.nlm.nih.gov/31367993/) → the stepwise treatment algorithm
130. Kremer AE, 2017, *Hautarzt* — Intrahepatic cholestasis of pregnancy: rare but important. [PMID 28074213](https://pubmed.ncbi.nlm.nih.gov/28074213/) → the clinical features of ICP itch (nocturnal · palms and soles)
131. Ovadia C, 2018, *J Clin Apher* — Therapeutic plasma exchange as a novel treatment for severe ICP. [PMID 30321466](https://pubmed.ncbi.nlm.nih.gov/30321466/) → the observation that plasma exchange works; the `NBD` node on the map (removing the autotaxin drive itself)

---

## 9. UDCA — trial results and mechanism

132. **Chappell LC, 2019, *Lancet*** — Ursodeoxycholic acid versus placebo in women with ICP (PITCHES): a randomised controlled trial. [PMID 31378395](https://pubmed.ncbi.nlm.nih.gov/31378395/) → **the trial reproduced in section E of the model.** Composite endpoint 24.8% versus 27.8%, itch −0.7 cm/10 cm
133. Chappell LC, 2018, *Trials* — PITCHES protocol. [PMID 30482254](https://pubmed.ncbi.nlm.nih.gov/30482254/) → the exact definition of the composite endpoint (the model's decomposition follows this definition)
134. Fleminger J, 2021, *BJOG* — UDCA in ICP: a secondary analysis of the PITCHES trial. [PMID 33063439](https://pubmed.ncbi.nlm.nih.gov/33063439/) → subgroup analysis; effect by severity
135. **Ovadia C, 2021, *Lancet Gastroenterol Hepatol*** — UDCA in ICP: systematic review and individual participant data meta-analysis. [PMID 33915090](https://pubmed.ncbi.nlm.nih.gov/33915090/) → **a result pointing in a different direction from PITCHES.** Section E of the model argues that the two studies merely measured different things and are not in contradiction, and sets out the arithmetic
136. Walker KF, 2020, *Cochrane Database Syst Rev* — Pharmacological interventions for treating ICP. [PMID 32716060](https://pubmed.ncbi.nlm.nih.gov/32716060/) → a systematic synthesis of effect sizes by drug
137. Kong X, 2016, *Medicine (Baltimore)* — Effectiveness and safety of UDCA in ICP: meta-analysis. [PMID 27749550](https://pubmed.ncbi.nlm.nih.gov/27749550/) → **the fall in total bile acids on UDCA (the comparator for the model's −34~−46%)**
138. Grand'Maison S, 2014, *J Obstet Gynaecol Can* — UDCA treatment for ICP: meta-analysis. [PMID 25184983](https://pubmed.ncbi.nlm.nih.gov/25184983/) → same purpose
139. Shen Y, 2019, *Clin Drug Investig* — Is it necessary to perform pharmacological interventions for ICP? Bayesian network meta-analysis. [PMID 30357607](https://pubmed.ncbi.nlm.nih.gov/30357607/) → the ranking between drugs (the comparator for section F of the model)
140. Zhang L, 2015, *Eur Rev Med Pharmacol Sci* — UDCA and SAMe in ICP: multi-centred randomized controlled trial. [PMID 26502869](https://pubmed.ncbi.nlm.nih.gov/26502869/) → the effect-size evidence for **`ESAM`, `ESAM_R` (SAMe)**
141. Zhang Y, 2016, *Hepat Mon* — UDCA and SAMe for ICP: meta-analysis. [PMID 27799965](https://pubmed.ncbi.nlm.nih.gov/27799965/) → same purpose
142. Marschall HU, 2015, *Expert Rev Gastroenterol Hepatol* — Management of ICP. [PMID 26313609](https://pubmed.ncbi.nlm.nih.gov/26313609/) → UDCA dose (10–20 mg/kg/d) → the model's dosing scenarios
143. Setoguchi T, 1984, *J Lipid Res* — Epimerization of the four 3,7-dihydroxy bile acid epimers by human fecal microorganisms. [PMID 6520544](https://pubmed.ncbi.nlm.nih.gov/6520544/) → the evidence for **`KDH_UD` (UDCA → LCA conversion)**
144. Hirano S, 1981, *J Lipid Res* — In vitro transformation of chenodeoxycholic acid and ursodeoxycholic acid by human intestinal flora. [PMID 7288282](https://pubmed.ncbi.nlm.nih.gov/7288282/) → same purpose, the magnitude of the conversion rate
145. Edenharder R, 1981, *J Lipid Res* — Epimerization of chenodeoxycholic acid to ursodeoxycholic acid by human intestinal Clostridia. [PMID 7276738](https://pubmed.ncbi.nlm.nih.gov/7276738/) → the reverse epimerisation (not included in the model — **an explicit simplification**)
146. Lee JY, 2013, *J Lipid Res* — Contribution of the 7β-hydroxysteroid dehydrogenase from *Ruminococcus gnavus* to UDCA formation in the human colon. [PMID 23729502](https://pubmed.ncbi.nlm.nih.gov/23729502/) → the organisms involved
147. Ferrari A, 1988, *Proc Soc Exp Biol Med* — In vitro transformation of cheno- and ursodeoxycholic acids by human intestinal microflora. [PMID 3368473](https://pubmed.ncbi.nlm.nih.gov/3368473/) → same purpose

---

## 10. Rifampicin · bile acid sequestrants · IBAT inhibitors · others

148. Hague WM, 2021, *BMC Pregnancy Childbirth* — A multi-centre, open label, randomised trial comparing UDCA with RIFampicin in ICP (TURRIFIC). [PMID 33435904](https://pubmed.ncbi.nlm.nih.gov/33435904/) → **the protocol of the trial that actually performs the comparison section F of the model predicts.** The model's prediction (that rifampicin beats UDCA substantially on itch and falls behind it on bile acids) is directly testable by this trial
149. Markus C, 2021, *Clin Chem Lab Med* — The BACH project: international total Bile Acid Comparison and Harmonisation project, sub-study of TURRIFIC. [PMID 34355544](https://pubmed.ncbi.nlm.nih.gov/34355544/) → that **the total bile acid assay is not harmonised between laboratories**; the practical reason the model reports the "assay value" and the "endogenous fraction" separately
150. Kremer AE, 2014, *Dig Dis* — Advances in pathogenesis and management of pruritus in cholestasis. [PMID 25034299](https://pubmed.ncbi.nlm.nih.gov/25034299/) → the place of rifampicin · colestyramine · naltrexone · sertraline in the treatment ladder
151. Mittal A, 2016, *Curr Probl Dermatol* — Cholestatic Itch Management. [PMID 27578083](https://pubmed.ncbi.nlm.nih.gov/27578083/) → same purpose
152. Peverelle M, 2026, *Aliment Pharmacol Ther* — Review: Ileal Bile Acid Transport (IBAT) inhibitors as an emerging treatment for cholestatic liver disease. [PMID 41953994](https://pubmed.ncbi.nlm.nih.gov/41953994/) → the pharmacological basis of the **`KIODEV` (IBAT inhibition) scenario**. There is no indication in pregnancy, and the model too marks it as hypothesis only
153. Lacey G, 2025, *Clin Ther* — Indirect comparison of maralixibat and odevixibat for PFIC. [PMID 40544071](https://pubmed.ncbi.nlm.nih.gov/40544071/) → the magnitude of the effect size
154. Muntaha HST, 2022, *J Clin Med* — IBAT blockers for cholestatic liver disease in Alagille syndrome: systematic review and meta-analysis. [PMID 36556142](https://pubmed.ncbi.nlm.nih.gov/36556142/) → the size of the fall in bile acids
155. Jeyaraj R, 2024, *Lancet Child Adolesc Health* — Paediatric research sets new standards for therapy in paediatric and adult cholestasis. [PMID 38006895](https://pubmed.ncbi.nlm.nih.gov/38006895/) → the place of the IBAT inhibitors
156. ElSalem SA, 2026, *Cureus* — Refractory ICP in twin gestation managed with UDCA and adjunctive cholestyramine. [PMID 41694888](https://pubmed.ncbi.nlm.nih.gov/41694888/) → a clinical case of adjunctive colestyramine
157. Pongcharoen P, 2016, *Eur J Pain* — An evidence-based review of systemic treatments for itch. [PMID 26416344](https://pubmed.ncbi.nlm.nih.gov/26416344/) → the level of evidence for naltrexone · sertraline · antihistamines → why `ENTX` and `EAH` are set small
158. Anderson SL, 2026, *Drugs Context* — New and emerging treatments for PBC-related pruritus. [PMID 42232646](https://pubmed.ncbi.nlm.nih.gov/42232646/) → the pipeline of new drugs on the itch axis
159. Levy C, 2023, *Clin Gastroenterol Hepatol* — New Treatment Paradigms in Primary Biliary Cholangitis. [PMID 36809835](https://pubmed.ncbi.nlm.nih.gov/36809835/) → the two-axis structure of cholestasis treatment (confirmation in another disease)
160. Dong XR, 2023, *World J Clin Cases* — Effect of polyene phosphatidylcholine/UDCA/ademetionine on pregnancy outcomes in ICP. [PMID 37900240](https://pubmed.ncbi.nlm.nih.gov/37900240/) → combination therapy
161. Liao E, 2025, *J Mol Histol* — Inhibition of the GPR30-PI3K pathway by 4-phenylbutyric acid in the treatment of ICP. [PMID 40063113](https://pubmed.ncbi.nlm.nih.gov/40063113/) → a new target (not included in the model)
162. Zeng Z, 2025, *J Reprod Immunol* — Probiotic VSL#3 alleviates ICP by upregulating FXR-FGF15. [PMID 41066868](https://pubmed.ncbi.nlm.nih.gov/41066868/) → intervention on the gut microbiota (an approach that takes the model's `KDH_*` and `EC50FXI` as its target)

---

## 11. Vitamin K · coagulation · postpartum haemorrhage

163. **Cemortan M, 2025, *BMC Pregnancy Childbirth*** — Comparative analysis of vitamin K levels in women with ICP. [PMID 40200185](https://pubmed.ncbi.nlm.nih.gov/40200185/) → the human data for **`MIC50`, `KVK`, `WINR` (vitamin K status and coagulation)**
164. Arslanoğlu T, 2025, *BMC Pregnancy Childbirth* — ICP and coagulation: a dual risk of hypercoagulability and bleeding. [PMID 40281473](https://pubmed.ncbi.nlm.nih.gov/40281473/) → the observation that **the coagulation abnormality of ICP runs in both directions**; the model represents only the bleeding side — **an explicit omission**
165. Mao J, 2025, *J Glob Health* — Association between perinatal complications and venous thromboembolism in postpartum women. [PMID 40375726](https://pubmed.ncbi.nlm.nih.gov/40375726/) → the thrombotic side of the risk (not in the model)

---

## 12. Twin pregnancy · recurrence · long-term prognosis

166. **Axelsen SM, 2024, *Acta Obstet Gynecol Scand*** — The effect of twin pregnancy in ICP: a case control study. [PMID 39058263](https://pubmed.ncbi.nlm.nih.gov/39058263/) → the evidence for **`TWIN = 1.55`**
167. Gu L, 2026, *J Matern Fetal Neonatal Med* — Maternal and neonatal outcomes of ICP in twin pregnancies: systematic review and meta-analysis. [PMID 42392869](https://pubmed.ncbi.nlm.nih.gov/42392869/) → outcomes of ICP in twins
168. Xu T, 2022, *BMC Pregnancy Childbirth* — Perinatal outcomes associated with ICP in twin pregnancies were worse than singletons. [PMID 36335293](https://pubmed.ncbi.nlm.nih.gov/36335293/) → same purpose
169. Zhao Y, 2025, *BMC Pregnancy Childbirth* — Preterm birth and stillbirth: total bile acid levels in ICP and outcomes of twin pregnancies. [PMID 40389846](https://pubmed.ncbi.nlm.nih.gov/40389846/) → stratification in twins
170. Mitta K, 2023, *Case Rep Womens Health* — Selective feticide reverses ICP in twins discordant for growth. [PMID 37534193](https://pubmed.ncbi.nlm.nih.gov/37534193/) → the observation that **reducing the placental load reverses the disease**; the most direct evidence supporting the structure in which `TWIN` acts only on the sex hormone trajectory
171. Rosenberg HM, 2026, *Obstet Gynecol* — ICP Recurrence in a Subsequent Pregnancy. [PMID 40811826](https://pubmed.ncbi.nlm.nih.gov/40811826/) → the recurrence rate; in the model, recurrence is a consequence of the genetic susceptibility vector being fixed
172. Zloto K, 2026, *Int J Gynaecol Obstet* — Outcomes of subsequent pregnancy following ICP. [PMID 42460752](https://pubmed.ncbi.nlm.nih.gov/42460752/) → same purpose
173. Sarker MR, 2024, *Am J Perinatol* — History of cholestasis is not associated with worsening outcomes in subsequent pregnancy. [PMID 38423120](https://pubmed.ncbi.nlm.nih.gov/38423120/) → severity does not necessarily worsen on recurrence
174. **Marschall HU, 2013, *Hepatology*** — ICP and associated hepatobiliary disease: a population-based cohort study. [PMID 23564560](https://pubmed.ncbi.nlm.nih.gov/23564560/) → the long-term risk of hepatobiliary disease (`E_LATER` on the map)
175. Wikström Shemer EA, 2015, *J Hepatol* — ICP and cancer, immune-mediated and cardiovascular diseases: population-based cohort study. [PMID 25772037](https://pubmed.ncbi.nlm.nih.gov/25772037/) → the breadth of the long-term prognosis
176. Odutola PO, 2023, *ILIVER* — ICP is associated with increased risk of hepatobiliary disease and adverse fetal outcomes: systematic review and meta-analysis. [PMID 40636920](https://pubmed.ncbi.nlm.nih.gov/40636920/) → same purpose
177. Wei Y, 2025, *Clin Rheumatol* — Postpartum may be a risk factor for biochemical flares in PBC. [PMID 40670882](https://pubmed.ncbi.nlm.nih.gov/40670882/) → the behaviour of cholestatic disease after delivery
178. Dai S, 2026, *Arch Gynecol Obstet* — ICP and offspring neurodevelopment: bile acid-mediated mechanisms and long-term outcomes. [PMID 41949636](https://pubmed.ncbi.nlm.nih.gov/41949636/) → long-term neurodevelopment of the offspring (not in the model — **a candidate for extension**)
179. Guner Yilmaz B, 2025, *Children (Basel)* — Early metabolic profile in neonates with maternal ICP. [PMID 41462795](https://pubmed.ncbi.nlm.nih.gov/41462795/) → neonatal metabolism

---

## 13. Comorbidities · modifying factors

180. Han Z, 2026, *Hepatol Commun* — Clinical outcomes of ICP with versus without chronic hepatitis B. [PMID 41493830](https://pubmed.ncbi.nlm.nih.gov/41493830/) → concomitant HBV (not in the model)
181. Gao Q, 2024, *BMC Pregnancy Childbirth* — ICP combined with different stages of HBV infection on pregnancy outcomes. [PMID 38582906](https://pubmed.ncbi.nlm.nih.gov/38582906/) → same purpose
182. Li X, 2024, *Diabetol Metab Syndr* — Gestational diabetes aggravates adverse perinatal outcomes in women with ICP. [PMID 38429774](https://pubmed.ncbi.nlm.nih.gov/38429774/) → concomitant GDM
183. Yang C, 2026, *Medicine (Baltimore)* — Association of GDM with ICP: Mendelian randomization. [PMID 41861211](https://pubmed.ncbi.nlm.nih.gov/41861211/) → the direction of causation
184. Lv M, 2026, *J Glob Health* — Impact of hypothyroidism on the risk of ICP. [PMID 41524244](https://pubmed.ncbi.nlm.nih.gov/41524244/) → the relationship with thyroid function
185. Kwon Y, 2025, *J Pharm Pharm Sci* — ICP associated with azathioprine: disproportionality analysis. [PMID 41477595](https://pubmed.ncbi.nlm.nih.gov/41477595/) → drug-induced ICP (the situation of lowering `GBSEP` pharmacologically)
186. Misirlioglu R, 2026, *Medicine (Baltimore)* — Prevalence and clinical significance of proteinuria in ICP. [PMID 42065144](https://pubmed.ncbi.nlm.nih.gov/42065144/) → concomitant proteinuria
187. Watad H, 2024, *PLoS One* — Proteinuria is a clinical characteristic of ICP but not a marker of severity. [PMID 39259746](https://pubmed.ncbi.nlm.nih.gov/39259746/) → same topic, unsuitable as a severity marker

---

## 14. Diagnostic prediction models (the model's comparators)

188. Han S, 2026, *J Matern Fetal Neonatal Med* — Innovative nomogram integrating bile acid metabolomics for early diagnosis of ICP. [PMID 42161851](https://pubmed.ncbi.nlm.nih.gov/42161851/) → **bile acid metabolomics (measurement by species) improves diagnosis** → clinical confirmation of why this QSP model computes species by species
189. Ding F, 2026, *BMC Pregnancy Childbirth* — Nomogram integrating dynamic total bile acid monitoring for preterm birth prediction in twin pregnancies. [PMID 42104311](https://pubmed.ncbi.nlm.nih.gov/42104311/) → **serial measurement is better than a single measurement** → the evidence behind the limitation section G of the model points out about itself (the smooth trajectory)
190. Xie W, 2025, *BMC Pregnancy Childbirth* — Nomogram for predicting preterm birth in ICP. [PMID 39984873](https://pubmed.ncbi.nlm.nih.gov/39984873/) → prediction of preterm birth
191. Wen M, 2025, *Medicine (Baltimore)* — Clinical prediction models for preterm birth in ICP. [PMID 41465959](https://pubmed.ncbi.nlm.nih.gov/41465959/) → same purpose
192. Yılmaz EBS, 2025, *BMC Pregnancy Childbirth* — Diagnostic and prognostic value of APRI and de Ritis ratio in ICP. [PMID 40877816](https://pubmed.ncbi.nlm.nih.gov/40877816/) → liver enzyme ratio indices
193. Kahraman NÇ, 2026, *J Perinat Med* — FAR, PAR, APRI and adverse neonatal outcomes in ICP. [PMID 41351215](https://pubmed.ncbi.nlm.nih.gov/41351215/) → same purpose
194. Gregorc P, 2025, *Diagnostics (Basel)* — Evaluating the effect of bile acid levels on maternal and perinatal outcomes in ICP. [PMID 40941672](https://pubmed.ncbi.nlm.nih.gov/40941672/) → external validation of the bile acid-outcome relationship

---

## 15. QSP methodology · mrgsolve

195. **Elmokadem A, 2019, *CPT Pharmacometrics Syst Pharmacol*** — Quantitative Systems Pharmacology and Physiologically-Based Pharmacokinetic Modeling With mrgsolve: A Hands-On Tutorial. [PMID 31652028](https://pubmed.ncbi.nlm.nih.gov/31652028/) → the implementation conventions this repository follows
196. Lu T, 2024, *CPT Pharmacometrics Syst Pharmacol* — gPKPDviz: A flexible R shiny tool for PK/PD simulations using mrgsolve. [PMID 38082557](https://pubmed.ncbi.nlm.nih.gov/38082557/) → the reference for the structure of the Shiny app
197. Nigam SK, 2020, *Clin Pharmacol Ther* — The Systems Biology of Drug Metabolizing Enzymes and Transporters: Relevance to QSP. [PMID 32119114](https://pubmed.ncbi.nlm.nih.gov/32119114/) → the methodological basis for transporter-centred QSP
198. Generaux G, 2019, *Pharmacol Res Perspect* — QST reproduces species differences in PF-04895162 liver safety due to combined mitochondrial and bile acid toxicity. [PMID 31624633](https://pubmed.ncbi.nlm.nih.gov/31624633/) → a precedent for treating bile acid toxicity with QSP
199. Woodhead JL, 2024, *Xenobiotica* — Prediction of the liver safety profile of a first-in-class myeloperoxidase inhibitor using QST modeling. [PMID 38874513](https://pubmed.ncbi.nlm.nih.gov/38874513/) → the same line of work
200. Mayo AK, 2026, *Clin Pharmacol Ther* — QST Model Predicts Obeticholic Acid-Associated Liver Injury in MASLD. [PMID 42332345](https://pubmed.ncbi.nlm.nih.gov/42332345/) → prediction of liver injury from a bile acid receptor agonist
201. Beaudoin JJ, 2023, *Int J Mol Sci* — Combination of a human biomimetic liver microphysiology system with a QST modeling platform. [PMID 37298645](https://pubmed.ncbi.nlm.nih.gov/37298645/) → coupling in vitro work with the model
202. Nyholm I, 2025, *J Hepatol* — Accumulation of altered serum bile acids predicts liver injury after portoenterostomy in biliary atresia. [PMID 39889904](https://pubmed.ncbi.nlm.nih.gov/39889904/) → confirmation in another disease that **bile acid composition (by species) predicts better than the total does**

---

## Where the model departs from the literature, or was set without it (the honest list)

This is the section a reviewer should read first.

| Item | Status |
|------|------|
| **Delivery timing in the 40–99 µmol/L band** | The model computes 39–40 weeks as optimal and **does not reproduce** the guidelines' recommendation of 37–38 weeks. Section G of `icp_calibration.py` works out the reason (the excess risk over this band in Ovadia is only 0.15%p) and the neonatal morbidity weight (<0.01) that would be needed for 37 weeks to become optimal. This may be a defect of the model, or it may mean that the basis of the recommendation does not lie in the stillbirth numbers; which of the two it is cannot be settled by the model alone |
| **`VP` (maximum placental transport capacity)** | No directly measured human value was found. A **model choice** fitted to the cord blood concentration of normal pregnancy |
| **`WTOX` cytotoxicity weights** | The **ordering** is supported by the hydrophobicity sequence and the cardiomyocyte experiments (references 93, 95, 124), but the **spacing** is a model choice |
| **`HLP` = 4 (the LPA→itch Hill coefficient)** | Not the cooperativity of a single binding event but the apparent value of the LPAR–TRPV1–GRPR multi-step amplification. Not to be interpreted mechanistically |
| **Direct fetal cardiac protection by UDCA** | References 96 and 97 report that UDCA protects the fetal myocardium directly. The model represents this only as `WTOX[UDCA]=0.02` and adds no separate term, so it **may underestimate the fetal protection given by UDCA** |
| **Enterohepatic recycling of progesterone sulfates** | Reported by reference 52, but the model keeps it only as a plasma state variable (simplification) |
| **Reverse CDCA↔UDCA epimerisation** | Reported by reference 145, but the model contains only UDCA→LCA (simplification) |
| **Thrombotic risk in ICP** | References 164 and 165 report both bleeding and thrombotic risk. The model represents only bleeding (vitamin K) |
| **Long-term neurodevelopment of the offspring** | Reference 178. Outside the scope of the model |
| **MRGPRX4 genotype** | Reference 118. A likely contributor to individual variation in itch, but not included in the model |
| **Within-day and between-day variation in maternal bile acids** | The model runs a smooth trajectory. Reference 189 reports that serial measurement is better, and this is the most important direction in which to extend the model |
| **Total number of fitted constants** | `HSB0`, `HSBSC`, `HN` (stillbirth, fitted to the three values in reference 14) · `HPT0`, `HM0`, `HMSC` (background rates) · `VSCALE` (the VAS scale) · `WMORB` (the utility ratio of neonatal morbidity to stillbirth, a judgement value) — eight in all. Everything else is either a transport, binding or turnover constant, or a scale factor set so that normal pregnancy lands on normal values |

---

## Guide to the files

| File | Contents |
|------|------|
| `icp_qsp_model.dot` / `.svg` / `.png` | 20-cluster mechanistic map |
| `icp_reference_model.py` | Dependency-free Python RK4 reference implementation (78 ODEs) — every equation was run here first |
| `icp_calibration.py` | Fitting, ablation and trial-reproduction analyses (sections A–I) |
| `icp_calibration_output.txt` | Output of the script above |
| `icp_mrgsolve_model.R` | mrgsolve model (78 ODEs, 35 scenarios) |
| `icp_shiny_app.R` | 10-tab Shiny dashboard |
| `README.md` | A summary of the model's structure and results |
