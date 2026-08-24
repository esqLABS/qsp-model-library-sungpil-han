# Chronic Hepatitis D (HDV) — QSP model references

> **Citation convention.** The links below are of two kinds.
> - **`PMID xxxxxxxx`** — a record **looked up directly on PubMed and confirmed** in this session.
>   The author, year, journal, and title have all been verified.
> - **`PubMed search`** — rather than guessing at a PMID that could not be confirmed, a **search link** that
>   will find the paper is given instead. Conference material (D-LIVR and the like) and
>   entries whose text could not be checked fall here.
>
> Whether a parameter value came from the literature below, was back-solved from a steady-state
> condition, or is a **hypothesis** of the model is stated item by item in sections A1/A3/A6/A13 of `hdv_model_report.txt`
> and in the comments of `hdv_reference_model.py`.

---

## 1. The receptor and entry (the NTCP entry axis) — the target of bulevirtide

1. Yan H, Zhong G, Xu G, et al. **Sodium taurocholate cotransporting polypeptide is a functional receptor for human hepatitis B and D virus.** *eLife* 2012;1:e00049. — the paper identifying NTCP (SLC10A1) as the HDV/HBV receptor. The whole entry term of this model rests on it. [PMID 23150796](https://pubmed.ncbi.nlm.nih.gov/23150796/)
2. Urban S, Neumann-Haefelin C, Lampertico P. **Hepatitis D virus in 2021: virology, immunology and new treatment approaches for a difficult-to-treat disease.** *Gut* 2021;70:1782-1794. — the standard review, integrating entry, replication, assembly, and immunity. [PMID 34103404](https://pubmed.ncbi.nlm.nih.gov/34103404/)
3. Mentha N, Clément S, Negro F, Alfaiate D. **A review on hepatitis D: From virology to new therapies.** *J Adv Res* 2019;17:3-15. [PMID 31193285](https://pubmed.ncbi.nlm.nih.gov/31193285/)
4. Ni Y, Lempp FA, Mehrle S, et al. **Hepatitis B and D viruses exploit sodium taurocholate co-transporting polypeptide for species-specific entry into hepatocytes.** *Gastroenterology* 2014. — the preS1-mimicry mechanism of myrcludex B (the precursor of bulevirtide). [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Ni+Lempp+hepatitis+B+D+sodium+taurocholate+species-specific+entry+hepatocytes)

## 2. Intracellular replication, RNA editing (ADAR1), L-HDAg

5. Casey JL. **RNA editing in hepatitis delta virus.** *Curr Top Microbiol Immunol* 2006. — the mechanism by which A→I editing at the amber/W site generates L-HDAg from S-HDAg. [PMID 16903221](https://pubmed.ncbi.nlm.nih.gov/16903221/)
6. Casey JL. **Control of ADAR1 editing of hepatitis delta virus RNAs.** *Curr Top Microbiol Immunol* 2012. [PMID 21732238](https://pubmed.ncbi.nlm.nih.gov/21732238/)
7. Jayan GC, Casey JL. **Increased RNA editing and inhibition of hepatitis delta virus replication by high-level expression of ADAR1 and ADAR2.** *J Virol* 2002;76:3819-3827. — the model's `f_edit` ↑ → L-HDAg ↑ → suppression of replication pathway. [PMID 11907222](https://pubmed.ncbi.nlm.nih.gov/11907222/)
8. Jayan GC, Casey JL. **Inhibition of hepatitis delta virus RNA editing by short inhibitory RNA-mediated knockdown of ADAR1 but not ADAR2.** *J Virol* 2002;76:12399-12404. [PMID 12414985](https://pubmed.ncbi.nlm.nih.gov/12414985/)
9. Chen R, Linnstaedt SD, Casey JL. **RNA editing and its control in hepatitis delta virus replication.** *Viruses* 2010;2:131-146. [PMID 21994604](https://pubmed.ncbi.nlm.nih.gov/21994604/)
10. Polson AG, Ley HL 3rd, Bass BL, Casey JL. **Hepatitis delta virus RNA editing is highly specific for the amber/W site and is suppressed by hepatitis delta antigen.** *Mol Cell Biol* 1998;18:1919-1926. — the negative feedback in which HDAg itself suppresses editing. [PMID 9528763](https://pubmed.ncbi.nlm.nih.gov/9528763/)
11. Jayan GC, Casey JL. **Effects of conserved RNA secondary structures on hepatitis delta virus genotype I RNA editing, replication, and virus production.** *J Virol* 2005;79:11187-11193. [PMID 16103170](https://pubmed.ncbi.nlm.nih.gov/16103170/)
12. Linnstaedt SD, Kasprzak WK, Shapiro BA, Casey JL. **The fraction of RNA that folds into the correct branched secondary structure determines hepatitis delta virus type 3 RNA editing levels.** *RNA* 2009;15:1177-1187. — the difference in editing biology of genotype 3 (this model is based on HDV-1). [PMID 19383766](https://pubmed.ncbi.nlm.nih.gov/19383766/)
13. Linnstaedt SD, Kasprzak WK, Shapiro BA, Casey JL. **The role of a metastable RNA secondary structure in hepatitis delta virus genotype III RNA editing.** *RNA* 2006;12:1521-1533. [PMID 16790843](https://pubmed.ncbi.nlm.nih.gov/16790843/)

## 3. Farnesylation and virion assembly — the target of lonafarnib

14. Otto JC, Casey PJ. **The hepatitis delta virus large antigen is farnesylated both in vitro and in animal cells.** *J Biol Chem* 1996;271:4569-4572. — the first identification of L-HDAg CXXQ (Cys211) farnesylation. [PMID 8617711](https://pubmed.ncbi.nlm.nih.gov/8617711/)
15. Glenn JS, Marsters JC Jr, Greenberg HB. **Use of a prenylation inhibitor as a novel antiviral agent.** *J Virol* 1998;72:9303-9306. — proof of concept that prenylation inhibition is an antiviral strategy. [PMID 9765479](https://pubmed.ncbi.nlm.nih.gov/9765479/)
16. Bordier BB, Ohkanda J, Liu P, et al. **In vivo antiviral efficacy of prenylation inhibitors against hepatitis delta virus.** *J Clin Invest* 2003;112:407-414. — the basis for the model's `I_FT → assembly blockade` term. [PMID 12897208](https://pubmed.ncbi.nlm.nih.gov/12897208/)
17. Heller T, Liang TJ. **Denying the wolf access to sheep's clothing.** *J Clin Invest* 2003;112:319-321. — the commentary on the paper above. The source of the central metaphor of this model, "cut off the supply of envelope". [PMID 12897197](https://pubmed.ncbi.nlm.nih.gov/12897197/)
18. Verrier ER, Weiss A, Bach C, et al. **Loss of hepatitis D virus infectivity upon farnesyl transferase inhibitor treatment associates with increasing RNA editing rates revealed by a new RT-ddPCR method.** *Antiviral Res* 2022;198:105250. — the observation that an FTI raises the editing rate (the dashed arrow on the map). [PMID 35051490](https://pubmed.ncbi.nlm.nih.gov/35051490/)
19. Liang YJ, Teng W, Chen CL, et al. **Statin inhibits large hepatitis delta antigen-Smad3-twist-mediated epithelial-to-mesenchymal transition and hepatitis D virus secretion.** *J Biomed Sci* 2020;27:65. — the influence of the mevalonate pathway on the assembly flux. [PMID 32434501](https://pubmed.ncbi.nlm.nih.gov/32434501/)

## 4. Persistence through cell division and regeneration — the "floor" of the model

20. Giersch K, Bhadra OD, Volz T, et al. **Hepatitis delta virus persists during liver regeneration and is amplified through cell division both in vitro and in vivo.** *Gut* 2019;68:150-157. — **the most important structural basis of this model.** HDV passes through hepatocyte division, is handed on to the daughter cells, and is amplified. That is what creates the floor an entry inhibitor cannot get below. [PMID 29217749](https://pubmed.ncbi.nlm.nih.gov/29217749/)
21. Winer BY, Gaska JM, Ding Q, et al. **Preclinical assessment of antiviral combination therapy in a genetically humanized mouse model for hepatitis delta virus infection.** *Sci Transl Med* 2018;10:eaap9328. — the preclinical basis for combining entry inhibition with IFN. [PMID 29950446](https://pubmed.ncbi.nlm.nih.gov/29950446/)

## 5. Mathematical models of viral kinetics (the line this model descends from directly)

22. Guedj J, Rotman Y, Cotler SJ, et al. **Understanding early serum hepatitis D virus and hepatitis B surface antigen kinetics during pegylated interferon-alpha therapy via mathematical modeling.** *Hepatology* 2014;60:1902-1910. — the first- and second-phase kinetics of pegIFN and the coupling to HBsAg. The basis for the model separating `eps_prod` (first phase) from `delta_IFN` (second phase). [PMID 25098971](https://pubmed.ncbi.nlm.nih.gov/25098971/)
23. Goyal A, Murray JM. **Dynamics of in vivo hepatitis D virus infection.** *J Theor Biol* 2016;398:9-19. — an ODE model treating the HDV-infected cell pool and the dependence on HBsAg. [PMID 27012516](https://pubmed.ncbi.nlm.nih.gov/27012516/)
24. Canini L, Koh C, Cotler SJ, et al. **Pharmacokinetics and pharmacodynamics modeling of lonafarnib in patients with chronic hepatitis delta virus infection.** *Hepatol Commun* 2017;1:288-292. — lonafarnib PK/PD modelling. The comparator for the `IC50_ft` and the size of the assembly-blockade effect in this model. [PMID 29404459](https://pubmed.ncbi.nlm.nih.gov/29404459/)
25. Shekhtman L, Cotler SJ, Hershkovich L, et al. **Modelling hepatitis D virus RNA and HBsAg dynamics during nucleic acid polymer monotherapy suggest rapid turnover of HBsAg.** *Sci Rep* 2020;10:7837. — the estimate that the half-life of circulating HBsAg is short (the model's `c_s`). [PMID 32398799](https://pubmed.ncbi.nlm.nih.gov/32398799/)
26. Zakh R, Churkin A, Bietsch W, et al. **A Mathematical Model for early HBV and -HDV Kinetics during Anti-HDV Treatment.** *Mathematics (Basel)* 2021;9:3323. [PMID 35282153](https://pubmed.ncbi.nlm.nih.gov/35282153/)
27. Maya S, Hershkovich L, Rotman Y, et al. **Hepatitis delta virus RNA decline post-inoculation in human NTCP transgenic mice is biphasic.** *mBio* 2023;14:e0060023. — the biexponential decay structure. [PMID 37436080](https://pubmed.ncbi.nlm.nih.gov/37436080/)
28. Maya S, et al. **Hepatitis delta virus RNA decline post inoculation in human NTCP transgenic mice is biphasic.** *bioRxiv* 2023 (preprint). [PMID 36824865](https://pubmed.ncbi.nlm.nih.gov/36824865/)
29. Dahari H, Shudo E, Ribeiro RM, Perelson AS. **Modelling hepatitis C virus kinetics: the relationship between the infected cell loss rate and the final slope of viral decay.** *Antivir Ther* 2009;14:459-464. — the methodological basis for "the final slope = the loss rate of infected cells". This is why an entry inhibitor in this model has no first phase and only a second. [PMID 19474480](https://pubmed.ncbi.nlm.nih.gov/19474480/)

## 6. Bulevirtide — the trials and the PK/bile acids

30. Wedemeyer H, Aleman S, Brunetto MR, et al. **A Phase 3, Randomized Trial of Bulevirtide in Chronic Hepatitis D.** *N Engl J Med* 2023;389:22-32. — **MYR301.** Combined response at 48 weeks 45% on 2 mg / 48% on 10 mg / 2% untreated, HDV RNA response 71%/76%/4%, ALT normalisation 51%/56%/12%. Two of this model's anchors (A3, A6) and the main comparator for its predictions. [PMID 37345876](https://pubmed.ncbi.nlm.nih.gov/37345876/)
31. Asselah T, Chulanov V, Lampertico P, et al. **Bulevirtide Combined with Pegylated Interferon for Chronic Hepatitis D.** *N Engl J Med* 2024;390:1935-1946. — **MYR204.** The superiority of the combination in sustained response after the end of treatment. The comparator for the synergy prediction A8 of this model. [PMID 38842520](https://pubmed.ncbi.nlm.nih.gov/38842520/)
32. Asselah T, Wedemeyer H, Aleman S, et al. **Bulevirtide Monotherapy Is Safe and Well Tolerated in Chronic Hepatitis Delta: An Integrated Safety Analysis of Bulevirtide Clinical Trials at Week 48.** *Liver Int* 2025. — the dose dependence and the asymptomatic nature of the rise in total bile acids. The safety data the bile acid anchor of A1 refers to. [PMID 39648559](https://pubmed.ncbi.nlm.nih.gov/39648559/)
33. Lampertico P, Anolli MP, Roulot D, Wedemeyer H. **Antiviral therapy for chronic hepatitis delta: new insights from clinical trials and real-life studies.** *Gut* 2025. — the observation that in real-world cohorts the long-term response to bulevirtide keeps increasing with time. The basis for the model predicting a gentle continued fall rather than a plateau. [PMID 39663120](https://pubmed.ncbi.nlm.nih.gov/39663120/)
34. Sandmann L, et al. **Treatment response to bulevirtide is linked to amelioration of portal hypertension in patients with chronic hepatitis D.** *JHEP Rep* 2025. — the clinical basis for the response carrying through to portal pressure and platelets (the model's `PLT` and `Fib` axes). [PMID 41674894](https://pubmed.ncbi.nlm.nih.gov/41674894/)
35. Zhu V, et al. **Evaluation of the drug-drug interaction potential of the novel hepatitis B and D virus entry inhibitor bulevirtide at OATP1B in healthy volunteers.** *Front Pharmacol* 2023;14:1128547. — data showing that bulevirtide inhibits OATP1B only weakly. The basis for holding the accessory bile acid uptake route `r_oatp` as a separate term in the model. [PMID 37089922](https://pubmed.ncbi.nlm.nih.gov/37089922/)
36. Bogomolov P, Alexandrov A, Voronkova N, et al. **Treatment of chronic hepatitis D with the entry inhibitor myrcludex B: First results of a phase Ib/IIa study.** *J Hepatol* 2016;65:490-498. — the first proof of concept. [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Bogomolov+chronic+hepatitis+D+entry+inhibitor+myrcludex+B+phase+Ib%2FIIa)
37. Blank A, Meier K, Urban S, Haefeli WE, Blank A. **Clinical pharmacology of the HBV/HDV entry inhibitor myrcludex B / bulevirtide — pharmacokinetics and bile acid changes.** — the size of the rise in total bile acids at 2 mg against 10 mg. The two anchors of A1. [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Blank+myrcludex+B+bulevirtide+pharmacokinetics+bile+acids)
38. European Medicines Agency. **Hepcludex (bulevirtide) — EPAR product information.** — the approved regimen of 2 mg SC once daily, and more than dose-proportional exposure. [EMA](https://www.ema.europa.eu/en/medicines/human/EPAR/hepcludex)

## 7. Interferon alfa / lambda

39. Wedemeyer H, Yurdaydìn C, Dalekos GN, et al. **Peginterferon plus adefovir versus either drug alone for hepatitis delta.** *N Engl J Med* 2011;364:322-331. — **HIDIT-1.** 48 weeks of pegIFN, with about 26-28% HDV RNA negative 24 weeks after the end of treatment. Adding adefovir was of no use (= consistent with the model structure in which NUCs do not work against HDV). [PMID 21268724](https://pubmed.ncbi.nlm.nih.gov/21268724/)
40. Heidrich B, Yurdaydìn C, Kabaçam G, et al. **Late HDV RNA relapse after peginterferon alpha-based therapy of chronic hepatitis delta.** *Hepatology* 2014;60:87-97. — late relapse on long-term follow-up. The basis for the model's term in which the recovery from exhaustion gradually reverses. [PMID 24585488](https://pubmed.ncbi.nlm.nih.gov/24585488/)
41. Etzion O, Hamid S, Lurie Y, et al. **Treatment of chronic hepatitis D with peginterferon lambda—the phase 2 LIMT-1 clinical trial.** *Hepatology* 2023;77:2093-2103. — a 36% response (5/14) 24 weeks after the end of treatment at 180 µg, with a problem of hyperbilirubinaemia. The lambda arm of the model (a receptor confined to the hepatocyte, no haematological toxicity, hepatic toxicity present). [PMID 36800850](https://pubmed.ncbi.nlm.nih.gov/36800850/)
42. Grabowski J, Yurdaydìn C, Zachou K, et al. **Hepatitis D virus-specific cytokine responses in patients with chronic hepatitis delta before and during interferon alfa-treatment.** *Liver Int* 2011;31:1395-1405. — the change in HDV-specific immune response during IFN treatment. The model's recovery-from-exhaustion term. [PMID 21762356](https://pubmed.ncbi.nlm.nih.gov/21762356/)
43. Mederacke I, Yurdaydin C, Grosshennig A, et al. **Renal function during treatment with adefovir plus peginterferon alfa-2a vs either drug alone in hepatitis B/D co-infection.** *J Viral Hepat* 2012;19:387-395. [PMID 22571900](https://pubmed.ncbi.nlm.nih.gov/22571900/)
44. Mederacke I, Bremer B, Heidrich B, et al. **Anti-HDV immunoglobulin M testing in hepatitis delta revisited: correlations with disease activity and response to pegylated interferon-alpha2a treatment.** *Antivir Ther* 2012;17:305-312. — the model's anti-HD IgM activity index. [PMID 22293066](https://pubmed.ncbi.nlm.nih.gov/22293066/)
45. Wedemeyer H, Yurdaydin C, Hardtke S, et al. **Peginterferon alfa-2a plus tenofovir disoproxil fumarate for hepatitis D (HIDIT-II): a randomised, placebo-controlled, phase 2 trial.** *Lancet Infect Dis* 2019;19:275-286. — **HIDIT-2.** Extending to 96 weeks did not greatly improve the response after the end of treatment. [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=HIDIT-II+peginterferon+alfa-2a+tenofovir+hepatitis+D+randomised)

## 8. Lonafarnib / ritonavir

46. Yurdaydin C, Keskin O, Yurdcu E, et al. **A phase 2 dose-finding study of lonafarnib and ritonavir with or without interferon alpha for chronic delta hepatitis.** *Hepatology* 2022;75:1551-1565. — **the LOWR HDV-2/4 series.** The strategy of raising the exposure of a low dose of lonafarnib by ritonavir boosting and so reducing gastrointestinal toxicity. [PMID 34860418](https://pubmed.ncbi.nlm.nih.gov/34860418/)
47. Yurdaydin C, Idilman R, Kalkan Ç, et al. **Optimizing lonafarnib treatment for the management of chronic delta hepatitis: The LOWR HDV-1 study.** *Hepatology* 2018;67:1224-1236. — the study showing the dose-limiting toxicity to be gastrointestinal. [PMID 29152762](https://pubmed.ncbi.nlm.nih.gov/29152762/)
48. Keskin O, Yurdaydin C. **Emerging drugs for hepatitis D.** *Expert Opin Emerg Drugs* 2023;28:1-11. [PMID 37096555](https://pubmed.ncbi.nlm.nih.gov/37096555/)
49. Etzion O, Asselah T, Yurdaydin C, et al. **D-LIVR: a phase 3 study of lonafarnib boosted with ritonavir with and without peginterferon alfa in chronic hepatitis delta.** — the phase 3 trial reporting a 48-week combined response of about 10% on lonafarnib/ritonavir, about 19% on the triple combination, and about 2% on placebo. Since the material is largely conference-based, a search link is given. [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=D-LIVR+lonafarnib+ritonavir+peginterferon+hepatitis+delta+phase+3)

## 9. Epidemiology, natural history, clinical outcomes

50. Stockdale AJ, Kreuels B, Henrion MYR, et al. **The global prevalence of hepatitis D virus infection: Systematic review and meta-analysis.** *J Hepatol* 2020;73:523-532. — an estimate of anti-HDV positivity among HBsAg-positive people. [PMID 32335166](https://pubmed.ncbi.nlm.nih.gov/32335166/)
51. Chen HY, Shen DT, Ji DZ, et al. **Prevalence and burden of hepatitis D virus infection in the global population: a systematic review and meta-analysis.** *Gut* 2019;68:512-521. [PMID 30228220](https://pubmed.ncbi.nlm.nih.gov/30228220/)
52. Gish RG, Wong RJ, Di Tanna GL, et al. **Association of hepatitis delta virus with liver morbidity and mortality: A systematic literature review and meta-analysis.** *Hepatology* 2024;79:1129-1140. — the increase in the risk of cirrhosis, decompensation, HCC, and death when HDV is present. The comparator for A11 of the model. [PMID 37870278](https://pubmed.ncbi.nlm.nih.gov/37870278/)
53. Stockdale AJ, Chaponda M, Beloukas A, et al. **Prevalence of hepatitis D virus infection in sub-Saharan Africa: a systematic review and meta-analysis.** *Lancet Glob Health* 2017;5:e992-e1003. [PMID 28911765](https://pubmed.ncbi.nlm.nih.gov/28911765/)
54. Wong RJ, Kaufman HW, Niles JK, et al. **Estimating the prevalence of hepatitis delta virus infection among adults in the United States: A meta-analysis.** *Liver Int* 2024;44:1715-1723. [PMID 38563728](https://pubmed.ncbi.nlm.nih.gov/38563728/)
55. Shen DT, Han J, Ji DZ, et al. **Epidemiology estimates of hepatitis D in individuals co-infected with human immunodeficiency virus and hepatitis B virus, 2002-2018: A systematic review and meta-analysis.** *J Viral Hepat* 2021;28:1057-1067. [PMID 33877742](https://pubmed.ncbi.nlm.nih.gov/33877742/)
56. Negro F, Lok AS. **Hepatitis D: A Review.** *JAMA* 2023;330:2376-2387. — the standard clinical review (diagnosis, treatment, natural history). [PMID 37943548](https://pubmed.ncbi.nlm.nih.gov/37943548/)
57. Degasperi E, Anolli MP, Lampertico P, et al. **Hepatitis D Virus Infection: Pathophysiology, Epidemiology and Treatment. Report From the Third Delta Cure Meeting 2024.** *Liver Int* 2025. — a summary of the latest development pipeline (siRNA, NAP, HBsAg-targeted). [PMID 40540405](https://pubmed.ncbi.nlm.nih.gov/40540405/)
58. Rizzetto M. **Chronic Hepatitis D; at a Standstill?** *Dig Dis* 2016;34:303-307. — the perspective of the author who first described HDV. [PMID 27170382](https://pubmed.ncbi.nlm.nih.gov/27170382/)
59. Abbas Z, Abbas M. **Management of hepatitis delta: Need for novel therapeutic options.** *World J Gastroenterol* 2015;21:9461-9465. [PMID 26327754](https://pubmed.ncbi.nlm.nih.gov/26327754/)
60. Alfaiate D, Clément S, Gomes D, Goossens N, Negro F. **Chronic hepatitis D and hepatocellular carcinoma: A systematic review and meta-analysis of observational studies.** *J Hepatol* 2020;73:533-539. — the increased HCC risk of HDV. The basis for the model's `bF_hcc` and `bV_hcc`. [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Alfaiate+chronic+hepatitis+D+hepatocellular+carcinoma+systematic+review+meta-analysis)

## 10. HBsAg-targeted therapy (cutting off the envelope supply)

61. Bazinet M, Pântea V, Cebotarescu V, et al. **Safety and efficacy of REP 2139 and pegylated interferon alfa-2a for treatment-naive patients with chronic hepatitis B virus and hepatitis D virus co-infection (REP 301 and REP 301-LTF).** *Lancet Gastroenterol Hepatol* 2017. — the strategy of lowering HBsAg with a nucleic acid polymer (NAP). [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Bazinet+REP+2139+pegylated+interferon+hepatitis+B+D+co-infection+REP+301)
62. **The application to HDV of HBsAg-targeted agents such as elebsiran (VIR-2218) and bepirovirsen.** — the approach of lowering HBsAg with an siRNA or ASO and so cutting off the envelope supply to HDV. The `siRNA` arm of A9 in this model. [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=elebsiran+OR+VIR-2218+OR+bepirovirsen+hepatitis+delta+HBsAg)

---

## Where each number comes from

| Model item | Type of source | Basis |
|---|---|---|
| NTCP as the HDV receptor | literature | #1, #2 |
| `Kd_ntcp`, `r_oatp` | **fitted (2)** | the two anchors for the fold rise in total bile acids in #32 and #37 → the occupancy is then **derived** by back-solving |
| The more than dose-proportional PK of bulevirtide | literature | #38 (implemented as saturable target-mediated disposition) |
| `dth_Id_immune` | **fitted (1)** | the 71% virological response rate at 48 weeks on 2 mg in #30 MYR301 |
| `ALT_base` | **fitted (1)** | the 51% ALT normalisation at 48 weeks on 2 mg in #30 MYR301 |
| Division-mediated persistence (the floor) | literature (structure) | #20 — the quantitative value is back-solved from steady state |
| Locality of regeneration, `loc_renew` | back-solved | from the observed stability of HBsAg |
| IFN first phase `eps_prod` / second phase `delta_IFN` | literature (structure) | #22, #29 |
| `IC50_ft`, assembly blockade | literature | #14-16, #24 |
| Intracellular RNA accumulation (the lonafarnib paradox) | a model prediction | a structural consequence, in the same direction as #18 |
| ADAR1 induction by IFN → editing ↑ | literature | #5-#7 |
| Recovery from exhaustion produces the sustained response | **a model hypothesis** | in the same direction as #40 and #42, with no quantitative basis |
| `kappa_fresh` (the ALT-RNA dissociation) | **a model hypothesis** | the minimum structure needed to explain the observation in #30. The falsifying condition is in report section A13 |
| Fibrosis progression and regression | literature | #52, #34 |
| HCC risk | literature | #52, #60 |

> The items marked as a **hypothesis** of the model are not values measured directly in the literature, and
> a **falsifier** for each is stated in section A13 of `hdv_model_report.txt`.
