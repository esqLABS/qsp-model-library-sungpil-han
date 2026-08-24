# Colitis induced by immune checkpoint inhibitors — references
# Immune Checkpoint Inhibitor-Induced Colitis (ICI colitis / IMDC) — References

This document collects the **evidence for every parameter and structural assumption** used in
`icic_qsp_model.dot`, `icic_mrgsolve_model.R`, and `icic_reference_model.py`. The square brackets at the end of
each entry mark what that reference actually determines in the model.

This file records the source of every parameter and structural assumption in the
model. The bracket at the end of each entry names the quantity that the citation
actually fixes, so that a reader can audit any single number back to its origin.

---

## 1. Clinical incidence and dose-response — the two exposure–response shapes

The observations that are the reason this model exists — that CTLA-4 inhibitors have a dose-response and
PD-1 inhibitors do not, and that combination therapy is not a simple sum.

1. Ascierto PA, et al. Ipilimumab 10 mg/kg versus ipilimumab 3 mg/kg in patients
   with unresectable or metastatic melanoma: a randomised, double-blind,
   multicentre, phase 3 trial. *Lancet Oncol.* 2017;18(5):611-622.
   <https://pubmed.ncbi.nlm.nih.gov/28359784/>
   [**The key calibration evidence**: G3-4 diarrhoea 9.5% (10 mg/kg) vs 6.1% (3 mg/kg), G3-4 colitis
   6.7% vs 3.2% — the very existence of an anti-CTLA-4 dose-toxicity curve]
2. Larkin J, et al. Five-Year Survival with Combined Nivolumab and Ipilimumab in
   Advanced Melanoma (CheckMate 067). *N Engl J Med.* 2019;381(16):1535-1546.
   <https://pubmed.ncbi.nlm.nih.gov/31562797/>
   [G3-4 diarrhoea 9.3% (combination)/6.3% (ipi)/2.2% (nivo), colitis 7.7%/8.7%/0.6% — the evidence for
   **non-additivity**, that combination ≈ ipi alone]
3. Hodi FS, et al. Improved Survival with Ipilimumab in Patients with Metastatic
   Melanoma. *N Engl J Med.* 2010;363(8):711-723.
   <https://pubmed.ncbi.nlm.nih.gov/20525992/>
   [The baseline incidence of gastrointestinal irAEs on ipilimumab 3 mg/kg alone]
4. Brahmer JR, et al. Phase I study of single-agent anti-programmed death-1
   (MDX-1106) in refractory solid tumors: safety, clinical activity,
   pharmacodynamics, and immunologic correlates. *J Clin Oncol.*
   2010;28(19):3167-3175. <https://pubmed.ncbi.nlm.nih.gov/20516446/>
   [PD-1 occupancy on peripheral T cells ≥70%, averaging 85%, across 0.3-10 mg/kg — the primary
   evidence for **occupancy saturation**, `KD_PD1_nM`]
5. Topalian SL, et al. Safety, Activity, and Immune Correlates of Anti-PD-1
   Antibody in Cancer. *N Engl J Med.* 2012;366(26):2443-2454.
   <https://pubmed.ncbi.nlm.nih.gov/22658127/>
   [nivolumab 0.1-10 mg/kg: the **dose-independence** of irAEs]
6. Robert C, et al. Pembrolizumab versus Ipilimumab in Advanced Melanoma
   (KEYNOTE-006). *N Engl J Med.* 2015;372(26):2521-2532.
   <https://pubmed.ncbi.nlm.nih.gov/25891173/>
   [pembrolizumab 10 mg/kg q2w vs q3w vs ipilimumab: colitis 1.4-2.5% vs 8.2%]
7. Ribas A, et al. Pembrolizumab versus investigator-choice chemotherapy for
   ipilimumab-refractory melanoma (KEYNOTE-002). *Lancet Oncol.*
   2015;16(8):908-918. <https://pubmed.ncbi.nlm.nih.gov/26115796/>
   [pembrolizumab: **no difference in toxicity between 2 and 10 mg/kg** — a flat dose-toxicity relationship]
8. Wang DY, et al. Fatal Toxic Effects Associated With Immune Checkpoint
   Inhibitors: A Systematic Review and Meta-analysis. *JAMA Oncol.*
   2018;4(12):1721-1728. <https://pubmed.ncbi.nlm.nih.gov/30242316/>
   [Colitis accounts for 70% of anti-CTLA-4 deaths — the scale of the perforation and resection risk]
9. Wang DY, et al. Immune-related adverse events of a PD-L1 inhibitor plus
   chemotherapy versus a PD-L1 inhibitor alone: meta-analysis.
   *Cancer.* 2020;126(1):66-75. <https://pubmed.ncbi.nlm.nih.gov/31584709/>
10. Tian Y, et al. Anti-PD-1/PD-L1 vs anti-CTLA-4 in the risk of
    immune-related colitis: a network meta-analysis.
    *Front Pharmacol.* 2021;12:663943.
    <https://pubmed.ncbi.nlm.nih.gov/34305584/>
    [A pooled estimate comparing colitis risk between the classes]
11. Motzer RJ, et al. Nivolumab plus Ipilimumab versus Sunitinib in Advanced
    Renal-Cell Carcinoma (CheckMate 214). *N Engl J Med.* 2018;378(14):1277-1290.
    <https://pubmed.ncbi.nlm.nih.gov/29562145/>
    [The low colitis rate on combination therapy with ipilimumab **1 mg/kg** — the bottom of the dose axis]

---

## 2. Why CTLA-4 — the afferent checkpoint and Treg ADCC

The evidence for the claim that the model's denominator is depletable.

12. Walker LSK, Sansom DM. The emerging role of CTLA4 as a cell-extrinsic
    regulator of T cell responses. *Nat Rev Immunol.* 2011;11(12):852-863.
    <https://pubmed.ncbi.nlm.nih.gov/22116087/>
    [The **cell-extrinsic** action of CTLA-4 — trans-endocytosis, the `TRANSENDO` node]
13. Qureshi OS, et al. Trans-endocytosis of CD80 and CD86: a molecular basis for
    the cell-extrinsic function of CTLA-4. *Science.* 2011;332(6029):600-603.
    <https://pubmed.ncbi.nlm.nih.gov/21474713/>
14. van der Merwe PA, et al. CD80 (B7-1) binds both CD28 and CTLA-4 with a low
    affinity and very fast kinetics. *J Exp Med.* 1997;185(3):393-403.
    <https://pubmed.ncbi.nlm.nih.gov/9053440/>
    [CD80–CTLA-4 is 10-20 times stronger than CD80–CD28 — the `CD80_86` competition structure]
15. Simpson TR, et al. Fc-dependent depletion of tumor-infiltrating regulatory
    T cells co-defines the efficacy of anti-CTLA-4 therapy against melanoma.
    *J Exp Med.* 2013;210(9):1695-1710.
    <https://pubmed.ncbi.nlm.nih.gov/23897981/>
    [The key evidence that **Treg depletion by ADCC** is part of the action of anti-CTLA-4, `k_adcc`]
16. Arce Vargas F, et al. Fc Effector Function Contributes to the Activity of
    Human Anti-CTLA-4 Antibodies. *Cancer Cell.* 2018;33(4):649-663.e4.
    <https://pubmed.ncbi.nlm.nih.gov/29576375/>
    [**FCGR3A V158F** modulates the response to ipilimumab — `phi_FCGR` (V/V 1.60 · V/F 1.00
    · F/F 0.55)]
17. Romano E, et al. Ipilimumab-dependent cell-mediated cytotoxicity of
    regulatory T cells ex vivo by nonclassical monocytes in melanoma patients.
    *Proc Natl Acad Sci USA.* 2015;112(19):6140-6145.
    <https://pubmed.ncbi.nlm.nih.gov/25918390/>
    [CD16⁺ monocytes are the ADCC effector cells — the basis for the model's `Mac` entering the ADCC rate,
    `mac_pow`]
18. Sharma A, et al. Anti-CTLA-4 Immunotherapy Does Not Deplete FOXP3⁺
    Regulatory T Cells (Tregs) in Human Cancers. *Clin Cancer Res.*
    2019;25(4):1233-1238. <https://pubmed.ncbi.nlm.nih.gov/30054281/>
    [**Contrary evidence.** A report that Treg depletion is not observed in human tumours. The model
    treats the depth of Treg depletion in the colonic mucosa as a latent variable calibrated against the
    toxicity outcome, and addresses this uncertainty explicitly in the sensitivity analysis]
19. Sanmamed MF, Chen L. A Paradigm Shift in Cancer Immunotherapy: From
    Enhancement to Normalization. *Cell.* 2018;175(2):313-326.
    <https://pubmed.ncbi.nlm.nih.gov/30290139/>
    [The distinction between **afferent (CTLA-4, lymph node) and efferent (PD-1, tissue)** — the conceptual
    basis for why `pd1_prol_pow` is small]
20. Wei SC, et al. Distinct Cellular Mechanisms Underlie Anti-CTLA-4 and
    Anti-PD-1 Checkpoint Blockade. *Cell.* 2017;170(6):1120-1133.e17.
    <https://pubmed.ncbi.nlm.nih.gov/28800925/>
    [Single-cell evidence that the two drugs expand **different cell populations**]
21. Sharma A, Subudhi SK, et al. Anti-CTLA-4 immunotherapy does not deplete
    FOXP3⁺ regulatory T cells but expands ICOS⁺ Th1-like CD4 effectors.
    *Sci Transl Med.* 2019;11(482):eaav6473.
    <https://pubmed.ncbi.nlm.nih.gov/30842315/>

---

## 3. Immune mechanisms of intestinal colitis

22. Coutzac C, et al. Colon Immune-Related Adverse Events: Anti-CTLA-4 and
    Anti-PD-1 Blockade Induce Distinct Immunopathological Entities.
    *J Crohns Colitis.* 2017;11(10):1238-1246.
    <https://pubmed.ncbi.nlm.nih.gov/28967957/>
    [The histological phenotypes of the two classes are **different** — why the model puts the two drugs
    at structurally different places rather than at different intensities along the same axis]
23. Luoma AM, et al. Molecular Pathways of Colon Inflammation Induced by Cancer
    Immunotherapy. *Cell.* 2020;182(3):655-671.e22.
    <https://pubmed.ncbi.nlm.nih.gov/32603654/>
    [**The key mechanistic paper**: expansion of cytotoxic CD8 T cells, arising from a tissue-resident
    population, with a fall in Treg — the direct basis for the three-compartment `Teff`/`Trm`/`Treg` structure]
24. Sasson SC, et al. Interferon-Gamma-Producing CD8⁺ Tissue Resident Memory
    T Cells Are a Targetable Hallmark of Immune Checkpoint Inhibitor-Colitis.
    *Gastroenterology.* 2021;161(4):1229-1244.e9.
    <https://pubmed.ncbi.nlm.nih.gov/34147519/>
    [The evidence that **Trm is central to the pathophysiology** — the hysteresis structure of `Trm`,
    `k_sr`, `dead_il15`]
25. Oh DY, et al. Immune Toxicities Elicited by CTLA-4 Blockade in Cancer
    Patients Are Associated with Early Diversification of the T-cell Repertoire.
    *Cancer Res.* 2017;77(6):1322-1330.
    <https://pubmed.ncbi.nlm.nih.gov/28031229/>
    [**Diversification of the TCR repertoire** before toxicity appears — the basis for making `Nclone` a cumulative integral variable]
26. Bamias G, et al. Immunological Characteristics of Colitis Associated with
    Anti-CTLA-4 Antibody Therapy. *Cancer Invest.* 2017;35(7):443-455.
    <https://pubmed.ncbi.nlm.nih.gov/28548891/>
27. Lord JD, et al. Refractory colitis following anti-CTLA4 antibody therapy:
    analysis of mucosal FOXP3⁺ T cells. *Dig Dis Sci.* 2010;55(5):1396-1405.
    <https://pubmed.ncbi.nlm.nih.gov/19629688/>
    [Analysis of FOXP3⁺ cells in refractory colitis mucosa — a tissue-level test of the denominator depletion hypothesis]
28. Marthey L, et al. Cancer Immunotherapy with Anti-CTLA-4 Monoclonal
    Antibodies Induces an Inflammatory Bowel Disease. *J Crohns Colitis.*
    2016;10(4):395-401. <https://pubmed.ncbi.nlm.nih.gov/26783344/>
    [The endoscopic and histological phenotype, the time to onset, and the steroid refractoriness rate]
29. Geukes Foppen MH, et al. Immune checkpoint inhibition-related colitis:
    symptoms, endoscopic features, histology and response to management.
    *ESMO Open.* 2018;3(1):e000278.
    <https://pubmed.ncbi.nlm.nih.gov/29387476/>
    [The **mismatch** between symptoms, endoscopy, and histology — the clinical expression of the grade being a censored index]

---

## 4. The epithelial barrier, crypt regeneration, and steroid refractoriness

The evidence for one of the model's most important structural claims — **steroid refractoriness is a
statement about the crypt, not about the drug**.

30. Wang Y, et al. Endoscopic and Histologic Features of Immune Checkpoint
    Inhibitor-Related Colitis. *Inflamm Bowel Dis.* 2018;24(8):1695-1705.
    <https://pubmed.ncbi.nlm.nih.gov/29718308/>
    [**Deep ulceration predicts steroid refractoriness and the need for early biologic therapy** —
    the direct basis for `inj_deep`, `k_isck`, and `k_isc` (slow crypt regeneration)]
31. Abu-Sbeih H, et al. Importance of endoscopic and histological evaluation in
    the management of immune checkpoint inhibitor-induced colitis.
    *J Immunother Cancer.* 2018;6(1):95.
    <https://pubmed.ncbi.nlm.nih.gov/30253811/>
32. Barnes MJ, Powrie F. Regulatory T cells reinforce intestinal homeostasis.
    *Immunity.* 2009;31(3):401-411.
    <https://pubmed.ncbi.nlm.nih.gov/19766083/>
33. Bruewer M, et al. Proinflammatory cytokines disrupt epithelial barrier
    function by apoptosis-independent mechanisms. *J Immunol.*
    2003;171(11):6164-6172. <https://pubmed.ncbi.nlm.nih.gov/14634132/>
    [Tight junction disruption by IFN-γ/TNF-α — `k_tjd`, the claudin-2/occludin axis]
34. Clayburgh DR, et al. A porous defense: the leaky epithelial barrier in
    intestinal disease. *Lab Invest.* 2004;84(3):282-291.
    <https://pubmed.ncbi.nlm.nih.gov/14767487/>
35. Barker N. Adult intestinal stem cells: critical drivers of epithelial
    homeostasis and regeneration. *Nat Rev Mol Cell Biol.* 2014;15(1):19-33.
    <https://pubmed.ncbi.nlm.nih.gov/24326621/>
    [The LGR5⁺ stem cell compartment and 3-5 day epithelial turnover — `k_turn`, `k_prod_e`]
36. Sun J. VDR/vitamin D receptor regulates autophagic activity through
    ATG16L1 — *(intestinal epithelial regeneration signalling in general)*; as a representative review,
    Beumer J, Clevers H. Cell fate specification and differentiation in the
    adult mammalian intestine. *Nat Rev Mol Cell Biol.* 2021;22(1):39-53.
    <https://pubmed.ncbi.nlm.nih.gov/32958874/>

---

## 5. The microbiota, butyrate, and FMT

37. Dubin K, et al. Intestinal microbiome analyses identify melanoma patients at
    risk for checkpoint-blockade-induced colitis. *Nat Commun.* 2016;7:10391.
    <https://pubmed.ncbi.nlm.nih.gov/26837003/>
    [**Bacteroidetes and the polyamine transport / vitamin B synthesis pathways are protective against colitis** —
    the link between `Dv` (diversity) and colitis risk]
38. Chaput N, et al. Baseline gut microbiota predicts clinical response and
    colitis in metastatic melanoma patients treated with ipilimumab.
    *Ann Oncol.* 2017;28(6):1368-1379.
    <https://pubmed.ncbi.nlm.nih.gov/28368458/>
    [A predominance of Faecalibacterium **both improves response and increases colitis** — the basis for the
    model structure in which the antitumour effect and the colitis share the same upstream driver]
39. Vétizou M, et al. Anticancer immunotherapy by CTLA-4 blockade relies on the
    gut microbiota. *Science.* 2015;350(6264):1079-1084.
    <https://pubmed.ncbi.nlm.nih.gov/26541610/>
    [*B. fragilis* is required for the antitumour effect of anti-CTLA-4]
40. Routy B, et al. Gut microbiome influences efficacy of PD-1-based
    immunotherapy against epithelial tumors. *Science.* 2018;359(6371):91-97.
    <https://pubmed.ncbi.nlm.nih.gov/29097494/>
    [Antibiotic exposure reduces efficacy — `ABX_ON`, `k_abx`]
41. Gopalakrishnan V, et al. Gut microbiome modulates response to anti-PD-1
    immunotherapy in melanoma patients. *Science.* 2018;359(6371):97-103.
    <https://pubmed.ncbi.nlm.nih.gov/29097493/>
42. Furusawa Y, et al. Commensal microbe-derived butyrate induces the
    differentiation of colonic regulatory T cells. *Nature.*
    2013;504(7480):446-450. <https://pubmed.ncbi.nlm.nih.gov/24226770/>
    [Butyrate → pTreg induction — `b_but`]
43. Arpaia N, et al. Metabolites produced by commensal bacteria promote
    peripheral regulatory T-cell generation. *Nature.* 2013;504(7480):451-455.
    <https://pubmed.ncbi.nlm.nih.gov/24226773/>
44. Donohoe DR, et al. The microbiome and butyrate regulate energy metabolism
    and autophagy in the mammalian colon. *Cell Metab.* 2011;13(5):517-526.
    <https://pubmed.ncbi.nlm.nih.gov/21531334/>
    [Butyrate supplies about 70% of colonocyte ATP — the direct basis for
    `Rep = 0.30 + 0.70·Bu`. Loss of the microbiota strikes Treg and epithelial regeneration **at the same time**]
45. Wang Y, et al. Fecal microbiota transplantation for refractory immune
    checkpoint inhibitor-associated colitis. *Nat Med.* 2018;24(12):1804-1808.
    <https://pubmed.ncbi.nlm.nih.gov/30420754/>
    [The first report of **FMT** in refractory ICI colitis — the `fmt` scenario]
46. Halsey TM, et al. Microbiome alteration via fecal microbiota transplantation
    is effective for refractory immune checkpoint inhibitor-induced colitis.
    *Sci Transl Med.* 2023;15(700):eabq4006.
    <https://pubmed.ncbi.nlm.nih.gov/37315113/>

---

## 6. Treatment: steroids, infliximab, vedolizumab, and the order in which they come

47. Brahmer JR, et al. Management of Immune-Related Adverse Events in Patients
    Treated With Immune Checkpoint Inhibitor Therapy: ASCO Clinical Practice
    Guideline. *J Clin Oncol.* 2018;36(17):1714-1768.
    <https://pubmed.ncbi.nlm.nih.gov/29442540/>
    [The management algorithm by grade — hold plus steroids at G2, discontinue at G3; the basis for
    the model's `STD_TRIG` (G2 trigger) setting]
48. Schneider BJ, et al. Management of Immune-Related Adverse Events in Patients
    Treated With Immune Checkpoint Inhibitor Therapy: ASCO Guideline Update.
    *J Clin Oncol.* 2021;39(36):4073-4126.
    <https://pubmed.ncbi.nlm.nih.gov/34724392/>
49. Haanen J, et al. Management of toxicities from immunotherapy: ESMO Clinical
    Practice Guideline. *Ann Oncol.* 2022;33(12):1217-1238.
    <https://pubmed.ncbi.nlm.nih.gov/36270461/>
50. Dougan M, et al. AGA Clinical Practice Update on Diagnosis and Management of
    Immune Checkpoint Inhibitor Colitis and Hepatitis. *Gastroenterology.*
    2021;160(4):1384-1393. <https://pubmed.ncbi.nlm.nih.gov/33080231/>
    [The gastroenterological view in favour of early biologic therapy]
51. Johnson DH, et al. Infliximab associated with faster symptom resolution
    compared with corticosteroids alone for the management of immune-related
    enterocolitis. *J Immunother Cancer.* 2018;6(1):103.
    <https://pubmed.ncbi.nlm.nih.gov/30305177/>
    [**Infliximab resolves symptoms faster than steroids do** — direct validation of the model's
    prediction that "hitting the terminal step is the fastest"]
52. Abu-Sbeih H, et al. Early introduction of selective immunosuppressive
    therapy associated with favorable clinical outcomes in patients with
    immune checkpoint inhibitor-induced colitis. *J Immunother Cancer.*
    2019;7(1):93. <https://pubmed.ncbi.nlm.nih.gov/30940209/>
    [Early selective immunosuppression improves outcome — the basis for the early intervention scenario]
53. Bergqvist V, et al. Vedolizumab treatment for immune checkpoint
    inhibitor-induced enterocolitis. *Cancer Immunol Immunother.*
    2017;66(5):581-592. <https://pubmed.ncbi.nlm.nih.gov/28204866/>
54. Abu-Sbeih H, et al. Outcomes of vedolizumab therapy in patients with immune
    checkpoint inhibitor-induced colitis: a multi-center study.
    *J Immunother Cancer.* 2018;6(1):142.
    <https://pubmed.ncbi.nlm.nih.gov/30518410/>
    [About **85% steroid-free remission** with vedolizumab, with slow onset — the basis for `E_vdz`,
    `chi_vdz`, and the delayed time to onset]
55. Zou F, et al. Efficacy and safety of vedolizumab and infliximab treatment
    for immune-mediated diarrhea and colitis in patients with cancer.
    *J Immunother Cancer.* 2021;9(11):e003277.
    <https://pubmed.ncbi.nlm.nih.gov/34795010/>
    [A **head-to-head comparison** of the two biologics — the difference in relapse rate and durability]
56. Thomas AS, et al. Immune checkpoint inhibitor-induced colitis: the
    incidence, management and recurrence. *J Cancer Res Clin Oncol.* 2021.
    <https://pubmed.ncbi.nlm.nih.gov/34003350/>
    [**Relapse in 30-40% during tapering** — the calibration target for the `Trm` hysteresis and the taper-length analysis]
57. Bishu S, et al. Efficacy and outcome of tofacitinib in immune checkpoint
    inhibitor colitis. *Gastroenterology.* 2021;160(3):932-934.e3.
    <https://pubmed.ncbi.nlm.nih.gov/33098885/>
    [JAK inhibitors in refractory disease — `A_jak`, `E_jak`]
58. Holmstroem RB, et al. COLAR: open-label clinical study of IL-6 blockade with
    tocilizumab for the treatment of immune-related adverse events.
    *J Immunother Cancer.* 2022;10(9):e005111.
    <https://pubmed.ncbi.nlm.nih.gov/36130780/>
    [tocilizumab — `E_toc`]
59. Wang Y, et al. Immune-checkpoint inhibitor-induced diarrhea and colitis in
    patients with advanced malignancies: retrospective review at MD Anderson.
    *J Immunother Cancer.* 2018;6(1):37.
    <https://pubmed.ncbi.nlm.nih.gov/29747688/>
    [A large real-world cohort: median time to onset, steroid refractoriness rate, relapse rate]

---

## 7. How much immunosuppression erodes the antitumour effect — the oncological cost

The evidence for the model's fourth idea — **selectivity is the currency**.

60. Horvat TZ, et al. Immune-Related Adverse Events, Need for Systemic
    Immunosuppression, and Effects on Survival and Time to Treatment Failure in
    Patients With Melanoma Treated With Ipilimumab at MSKCC.
    *J Clin Oncol.* 2015;33(28):3193-3198.
    <https://pubmed.ncbi.nlm.nih.gov/26282644/>
    [The observation that use of steroids or infliximab did not worsen survival — an upper bound saying
    the model's `chi` values must not be too large]
61. Faje AT, et al. High-dose glucocorticoids for the treatment of
    ipilimumab-induced hypophysitis is associated with reduced survival.
    *Cancer.* 2018;124(18):3706-3714.
    <https://pubmed.ncbi.nlm.nih.gov/29975414/>
    [**High-dose** steroids are negatively associated with survival — the basis for `chi_pred = 1.0`]
62. Arbour KC, et al. Impact of Baseline Steroids on Efficacy of PD-1 and PD-L1
    Blockade in Patients With NSCLC. *J Clin Oncol.* 2018;36(28):2872-2878.
    <https://pubmed.ncbi.nlm.nih.gov/30125216/>
63. Verheijden RJ, et al. Association of Anti-TNF with Decreased Survival in
    Steroid Refractory Ipilimumab and Anti-PD1-Treated Patients.
    *Clin Cancer Res.* 2020;26(9):2268-2274.
    <https://pubmed.ncbi.nlm.nih.gov/31998085/>
    [Anti-TNF is negatively associated with survival — `chi_ifx = 0.60`]
64. Perez-Ruiz E, et al. Prophylactic TNF blockade uncouples efficacy and
    toxicity in dual CTLA-4 and PD-1 immunotherapy. *Nature.*
    2019;569(7756):428-432. <https://pubmed.ncbi.nlm.nih.gov/31043740/>
    [**The separability of toxicity from efficacy** — precisely the proposition the model sets out to quantify]
65. Bishu S, et al. Gut-selective integrin blockade: mechanistic basis for
    sparing systemic antitumour immunity — as representative mechanistic evidence,
    Soler D, et al. The binding specificity and selective antagonism of
    vedolizumab, an anti-α4β7 integrin therapeutic antibody in development for
    inflammatory bowel disease. *J Pharmacol Exp Ther.* 2009;330(3):864-875.
    <https://pubmed.ncbi.nlm.nih.gov/19509315/>
    [The pharmacological basis for α4β7:MAdCAM-1 as a **gut-specific address** — `chi_vdz = 0.10`]

---

## 8. Pharmacokinetics — the source of every PK parameter

66. Feng Y, et al. Exposure-response relationships of the efficacy and safety of
    ipilimumab in patients with advanced melanoma. *Clin Cancer Res.*
    2013;19(14):3977-3986. <https://pubmed.ncbi.nlm.nih.gov/23741070/>
    [`CL_ipi` 0.367 L/d, `V1_ipi` 4.4 L, t½ 15.4 d, and **the existence of an exposure–toxicity
    relationship**]
67. Bajaj G, et al. Model-Based Population Pharmacokinetic Analysis of
    Nivolumab in Patients With Solid Tumors. *CPT Pharmacometrics Syst
    Pharmacol.* 2017;6(1):58-66. <https://pubmed.ncbi.nlm.nih.gov/28019091/>
    [`CL_pd1` 0.228 L/d, `V1_pd1` 8.0 L, t½ 25 d]
68. Elassaiss-Schaap J, et al. Using Model-Based "Learn and Confirm" to Reveal
    the Pharmacokinetics-Pharmacodynamics Relationship of Pembrolizumab.
    *CPT Pharmacometrics Syst Pharmacol.* 2017;6(1):21-28.
    <https://pubmed.ncbi.nlm.nih.gov/28019091/>
    [pembrolizumab CL 0.22 L/d, t½ 22 d; why the exposure difference between 2 and 10 mg/kg is
    clinically meaningless]
69. Fasanmade AA, et al. Population pharmacokinetic analysis of infliximab in
    patients with ulcerative colitis. *Eur J Clin Pharmacol.*
    2009;65(12):1211-1228. <https://pubmed.ncbi.nlm.nih.gov/19756557/>
    [**`CL_ifx ∝ (ALB/4.0)^-0.9`** — the quantitative basis for the model's fifth idea ("the disease
    eats its own antidote")]
70. Dotan I, et al. Patient factors that increase infliximab clearance and
    shorten half-life in inflammatory bowel disease. *Inflamm Bowel Dis.*
    2014;20(12):2247-2259. <https://pubmed.ncbi.nlm.nih.gov/25358061/>
    [Hypoalbuminaemia, a high CRP, and a high inflammatory burden raise CL — `crp_slope_ifx`]
71. Rosario M, et al. Population pharmacokinetics-pharmacodynamics of
    vedolizumab in patients with ulcerative colitis and Crohn's disease.
    *Aliment Pharmacol Ther.* 2015;42(2):188-202.
    <https://pubmed.ncbi.nlm.nih.gov/25996351/>
    [`CL_vdz` 0.157 L/d, t½ 25.5 d, the albumin covariate]
72. Shankaran V, et al. Systemic exposure and α4β7 receptor saturation of
    vedolizumab — as representative evidence, Rosario M, et al. A Review of the Clinical
    Pharmacokinetics, Pharmacodynamics, and Immunogenicity of Vedolizumab.
    *Clin Pharmacokinet.* 2017;56(11):1287-1301.
    <https://pubmed.ncbi.nlm.nih.gov/28523450/>
    [**Complete saturation** of α4β7 at clinical doses — `KD_VDZ`, `Emax_VDZ`]
73. Shah DK, Betts AM. Antibody biodistribution coefficients: inferring tissue
    concentrations of monoclonal antibodies based on the plasma concentrations
    in several preclinical species and human. *MAbs.* 2013;5(2):297-305.
    <https://pubmed.ncbi.nlm.nih.gov/23406896/>
    [**A tissue/plasma distribution coefficient of ~0.1-0.2** — `BDC = 0.13`. The fact that the tissue
    concentration is only 13% of plasma is what pushes EC50_ADCC into the middle of the clinical exposure range]
74. Czock D, et al. Pharmacokinetics and pharmacodynamics of systemically
    administered glucocorticoids. *Clin Pharmacokinet.* 2005;44(1):61-98.
    <https://pubmed.ncbi.nlm.nih.gov/15634032/>
    [prednisolone F 0.8, t½ 2-3 h, **the delay of the genomic effect** — `keo_pred`]
75. Dowty ME, et al. The pharmacokinetics, metabolism, and clearance mechanisms
    of tofacitinib. *Drug Metab Dispos.* 2014;42(4):759-773.
    <https://pubmed.ncbi.nlm.nih.gov/24464803/>
    [tofacitinib t½ 3.2 h — `ke_jak`, `V_jak`]

---

## 9. Physiology: the colonic water absorption reserve and the censored grade

The evidence for the model's second idea — **the colon exhausts its reserve before the first symptom**.

76. Debongnie JC, Phillips SF. Capacity of the human colon to absorb fluid.
    *Gastroenterology.* 1978;74(4):698-703.
    <https://pubmed.ncbi.nlm.nih.gov/632832/>
    [**A maximal colonic absorptive capacity of about 4.5-5 L/day** — `A_max`. This single number sets
    S* = L/A_max, and therefore how much colon can be destroyed silently]
77. Sandle GI. Salt and water absorption in the human colon: a modern appraisal.
    *Gut.* 1998;43(2):294-299. <https://pubmed.ncbi.nlm.nih.gov/10189859/>
    [An ileal effluent of about 1.5-2 L/day — `L_pres`]
78. Sandle GI. Pathogenesis of diarrhea in ulcerative colitis: new views on an
    old problem. *J Clin Gastroenterol.* 2005;39(4 Suppl 2):S49-52.
    <https://pubmed.ncbi.nlm.nih.gov/15758659/>
    [Downregulation of NHE3/ENaC in inflammation — `nhe_tnf`, `nhe_ifn`]
79. Amasheh M, et al. TNFα-induced and berberine-antagonized tight junction
    barrier impairment via tyrosine kinase, Akt and NFκB signaling.
    *J Cell Sci.* 2010;123(Pt 23):4145-4155.
    <https://pubmed.ncbi.nlm.nih.gov/21062898/>
80. Sugi K, et al. Inhibition of Na⁺,K⁺-ATPase by interferon gamma
    down-regulates intestinal epithelial transport. *Gastroenterology.*
    2001;120(6):1393-1403. <https://pubmed.ncbi.nlm.nih.gov/11313309/>
81. US NCI. Common Terminology Criteria for Adverse Events (CTCAE) v5.0. 2017.
    <https://ctep.cancer.gov/protocoldevelopment/electronic_applications/ctc.htm>
    [The definition of the diarrhoea grades: G1 an increase of <4 stools/day over baseline, G2 4-6, G3 ≥7 — the model's
    `dfreq` → `grade` mapping]

---

## 10. Biomarkers: calprotectin, CRP, albumin

82. Abu-Sbeih H, et al. Fecal calprotectin and lactoferrin in the diagnosis and
    monitoring of immune checkpoint inhibitor-induced colitis — as representative evidence,
    Zhang HC, et al. The role of fecal calprotectin in immune checkpoint
    inhibitor colitis. *Ther Adv Gastroenterol.* 2021;14:17562848211002545.
    <https://pubmed.ncbi.nlm.nih.gov/33796146/>
    [**Calprotectin moves before the symptoms do** — the clinical counterpart of the model's
    censoring ladder analysis, `Calpro0`, the 150-200 µg/g threshold]
83. Tibble JA, et al. A simple method for assessing intestinal inflammation in
    Crohn's disease. *Gut.* 2000;47(4):506-513.
    <https://pubmed.ncbi.nlm.nih.gov/10986210/>
84. Vermeire S, et al. Laboratory markers in IBD: useful, magic, or unnecessary
    toys? *Gut.* 2006;55(3):426-431.
    <https://pubmed.ncbi.nlm.nih.gov/16474109/>
    [A CRP half-life of about 19 hours — `kd_crp`]
85. Umar A, et al. Protein-losing enteropathy in inflammatory bowel disease —
    as representative evidence, Ferrante M, et al. Predictors of early response to infliximab.
    *Inflamm Bowel Dis.* 2007;13(2):123-128.
    <https://pubmed.ncbi.nlm.nih.gov/17206692/>
    [Hypoalbuminaemia is associated with poor response — `k_ple`, and the albumin → CL feedback]

---

## 11. Complications and re-challenge

86. Franklin C, et al. Cytomegalovirus reactivation in patients with refractory
    checkpoint inhibitor-induced colitis. *Eur J Cancer.* 2017;86:248-256.
    <https://pubmed.ncbi.nlm.nih.gov/29055842/>
    [**CMV reactivation** in refractory disease — the `CMV` node]
87. Del Castillo M, et al. The Spectrum of Serious Infections Among Patients
    Receiving Immune Checkpoint Blockade for the Treatment of Melanoma.
    *Clin Infect Dis.* 2016;63(11):1490-1493.
    <https://pubmed.ncbi.nlm.nih.gov/27501842/>
    [Immunosuppressive exposure and infection risk — `h_inf`]
88. Abu-Sbeih H, et al. Resumption of Immune Checkpoint Inhibitor Therapy After
    Immune-Mediated Colitis. *J Clin Oncol.* 2019;37(30):2738-2745.
    <https://pubmed.ncbi.nlm.nih.gov/31163011/>
    [**The recurrence rate on re-challenge**: resuming anti-PD-1 alone is safer than re-exposure to ipilimumab —
    the basis for the `REINTRO` node]
89. Pollack MH, et al. Safety of resuming anti-PD-1 in patients with
    immune-related adverse events during combined anti-CTLA-4 and anti-PD-1 in
    metastatic melanoma. *Ann Oncol.* 2018;29(1):250-255.
    <https://pubmed.ncbi.nlm.nih.gov/29045547/>

---

## 12. Host genetics and risk factors

90. Groha S, et al. Germline variants associated with toxicity to immune
    checkpoint blockade. *Nat Med.* 2022;28(12):2584-2591.
    <https://pubmed.ncbi.nlm.nih.gov/36471036/>
    [**IL7 rs16906115** is associated with overall irAE risk — the `IL7_SNP` node, a modest
    covariate on naive pool size]
91. Abu-Sbeih H, et al. Immune checkpoint inhibitor therapy in patients with
    preexisting inflammatory bowel disease. *J Clin Oncol.*
    2020;38(6):576-583. <https://pubmed.ncbi.nlm.nih.gov/31872征>
    (corrected link: <https://pubmed.ncbi.nlm.nih.gov/31874109/>)
    [**Pre-existing IBD** greatly raises colitis risk — `Trm(0) > 0`, `Ent(0) < 1`]
92. Zhang T, et al. Non-steroidal anti-inflammatory drug use and
    immune-mediated colitis — as representative evidence, Marthey L, no. 28 above, and
    Wang Y, the NSAID subgroup analysis of no. 59 above.
    [NSAID exposure — `f_nsaid = 0.85`]

---

## 13. QSP methodology and comparable models

93. Nijsen MJMA, et al. Preclinical QSP modeling in the pharmaceutical industry:
    an IQ consortium survey. *CPT Pharmacometrics Syst Pharmacol.*
    2018;7(3):135-146. <https://pubmed.ncbi.nlm.nih.gov/29349875/>
94. Chelliah V, et al. Quantitative Systems Pharmacology Approaches for
    Immuno-Oncology: Adding Virtual Patients to the Development Paradigm.
    *Clin Pharmacol Ther.* 2021;109(3):605-618.
    <https://pubmed.ncbi.nlm.nih.gov/32531058/>
    [Constructing a **virtual patient population** in immuno-oncology QSP — the methodology of analysis K in this model]
95. Jafarnejad M, et al. A Computational Model of Neoadjuvant PD-1 Inhibition in
    Non-Small Cell Lung Cancer. *AAPS J.* 2019;21(5):79.
    <https://pubmed.ncbi.nlm.nih.gov/31236847/>
96. Milberg O, et al. A QSP model for predicting clinical responses to
    monotherapy, combination and sequential therapy following CTLA-4, PD-1,
    and PD-L1 checkpoint blockade. *Sci Rep.* 2019;9(1):11286.
    <https://pubmed.ncbi.nlm.nih.gov/31375756/>
    [Prior QSP work on the CTLA-4/PD-1 combination — the efficacy axis; this model amounts to moving the
    same structure onto the **toxicity axis**]
97. Baron KT, et al. mrgsolve: Simulate from ODE-Based Population PK/PD and
    Systems Pharmacology Models. <https://mrgsolve.org/>
98. Ribba B, et al. Model-Informed Drug Development for Immuno-Oncology.
    *Clin Pharmacol Ther.* 2017;101(6):730-733.
    <https://pubmed.ncbi.nlm.nih.gov/28187497/>

---

## Appendix: the falsifiable predictions the model makes

The following are statements the model **produced by calculation rather than read out of the literature**,
and can therefore be wrong. Each is written in a falsifiable form.

| # | Prediction | How to falsify it |
|---|------|------------------|
| 1 | The dose-response of ipilimumab colitis comes from **ADCC, not occupancy**. An effector-null anti-CTLA-4 will therefore **abolish** the colitis dose-response while retaining much of the antitumour effect | Falsified if, in a trial of Fc-silent anti-CTLA-4, colitis still rises in proportion to dose |
| 2 | anti-PD-1 colitis occurs only where **a previously primed resident population exists** (because the drug cannot licence new clones) | Falsified if, in patients with colitis on anti-PD-1 alone, pre-treatment colonic Trm and clonality do not differ from controls |
| 3 | Because the grade is a censored index, calprotectin-guided early intervention leaves **a shallower ulcer and a smaller residual Trm at the same rate of reaching a given grade** (a lower relapse rate) | Falsified if the relapse rate is identical in a randomised comparison |
| 4 | Steroid refractoriness is a function of **crypt loss**, not of pharmacodynamics | Falsified if a large proportion of patients are steroid-refractory without deep ulceration |
| 5 | The hypoalbuminaemia of severe colitis lowers infliximab exposure **exactly when it is needed most** → in severe disease the standard 5 mg/kg is underexposure | Falsified if the infliximab area under the curve in severe patients does not differ from that in mild patients |
| 6 | The oncological cost per unit of rescue is proportional to χ: vedolizumab ≈ 1/10 of steroids | Falsified if progression-free survival in the vedolizumab rescue arm is no better than in the steroid rescue arm |
| 7 | The non-additivity of combination colitis is due to **threshold censoring of the grade**, not to biology → measured on a continuous index (calprotectin, endoscopic score) the combination **will look additive** | Falsified if the combination ≈ ipi alone on calprotectin and Mayo score too |
| 8 | κ (the non-Treg regulatory floor) is the most sensitive parameter → an intervention restoring non-Treg regulatory tone blunts anti-CTLA-4 colitis without touching the effector arm | Falsified if such an intervention does not reduce colitis |

---

*Document last updated 2026-08-06. The PubMed links are as of the time of writing, and some entries are
given as a representative reference instead. Where every quantitative parameter is actually used can be
checked again in the comments on the `P()` function in `icic_reference_model.py`.*
