# Postpartum depression (PPD) QSP model — references
# References for the Postpartum Depression QSP Model

**Every PMID was looked up directly with NCBI E-utilities (`esearch` + `esummary`) in this session and its
title, journal, and year confirmed.** None of it was written from memory, and the verification script drew
on the same sources as the calibration targets recorded in the header of
`ppd_reference_check.py`. References that could not be confirmed were not put in this list.

Every PMID below was resolved live against PubMed during model construction; titles,
journals and years are as returned by NCBI. Nothing here is quoted from memory.

---

## 1. Epidemiology and clinical definition

1. Munk-Olsen T, et al. **New parents and mental disorders: a population-based register
   study.** *JAMA* 2006. — the temporal concentration of the risk of psychiatric admission after childbirth (the first 3 months).
   <https://pubmed.ncbi.nlm.nih.gov/17148723/>
2. Wisner KL, et al. **Onset timing, thoughts of self-harm, and diagnoses in postpartum
   women with screen-positive depression findings.** *JAMA Psychiatry* 2013. — a screening cohort of 10,000;
   the distribution of onset times and the frequency of self-harm ideation. The comparator for the model's "onset window".
   <https://pubmed.ncbi.nlm.nih.gov/23487258/>
3. O'Hara MW, McCabe JE. **Postpartum depression: current status and future
   directions.** *Annu Rev Clin Psychol* 2013.
   <https://pubmed.ncbi.nlm.nih.gov/23394227/>
4. Stewart DE, Vigod S. **Postpartum Depression.** *N Engl J Med* 2016.
   <https://pubmed.ncbi.nlm.nih.gov/27959754/>
5. Meltzer-Brody S, et al. **Postpartum psychiatric disorders.** *Nat Rev Dis Primers*
   2018. <https://pubmed.ncbi.nlm.nih.gov/29695824/>
6. Howard LM, et al. **Non-psychotic mental disorders in the perinatal period.**
   *Lancet* 2014. <https://pubmed.ncbi.nlm.nih.gov/25455248/>
7. Gavin NI, et al. **Perinatal depression: a systematic review of prevalence and
   incidence.** *Obstet Gynecol* 2005.
   <https://pubmed.ncbi.nlm.nih.gov/16260528/>
8. Shorey S, et al. **Prevalence and incidence of postpartum depression among healthy
   mothers: a systematic review and meta-analysis.** *J Psychiatr Res* 2018.
   <https://pubmed.ncbi.nlm.nih.gov/30114665/>
9. **Prevalence of Postpartum Depression Based on Diagnostic Interviews: A Systematic
   Review and Meta-Analysis.** *Depress Anxiety* 2023. — prevalence on structured interview criteria
   (lower than on screening criteria).
   <https://pubmed.ncbi.nlm.nih.gov/40224605/>
10. **Association of Postpartum Depression with Maternal Suicide: A Nationwide
    Population-Based Study.** *Int J Environ Res Public Health* 2022. — the justification for the model's
    suicidal ideation endpoint.
    <https://pubmed.ncbi.nlm.nih.gov/35564525/>
11. Bergink V, et al. **Postpartum Psychosis: Madness, Mania, and Melancholia in
    Motherhood.** *Am J Psychiatry* 2016. — a separate emergency condition, distinct from PPD.
    <https://pubmed.ncbi.nlm.nih.gov/27609245/>

## 2. Rating scales and their model mapping

12. Cox JL, Holden JM, Sagovsky R. **Detection of postnatal depression. Development of
    the 10-item Edinburgh Postnatal Depression Scale.** *Br J Psychiatry* 1987. — EPDS
    the original paper. <https://pubmed.ncbi.nlm.nih.gov/3651732/>
13. Levis B, et al. **Accuracy of the Edinburgh Postnatal Depression Scale (EPDS) for
    screening to detect major depression among pregnant and postpartum women.** *BMJ*
    2020. — an individual patient data meta-analysis; the basis for the 12/13 cut-off.
    <https://pubmed.ncbi.nlm.nih.gov/33177069/>
14. Zimmerman M, et al. **Is the cutoff to define remission on the Hamilton Rating
    Scale for Depression too high?** *J Nerv Ment Dis* 2005. — the basis for the model using HAM-D ≤ 7 as
    the boundary of remission, and the limits of doing so.
    <https://pubmed.ncbi.nlm.nih.gov/15729106/>
15. **Associations between commonly used patient-reported outcome tools in postpartum
    depression clinical practice and the Hamilton Rating Scale.** *Arch Womens Ment
    Health* 2020. — conversion between EPDS and HAM-D.
    <https://pubmed.ncbi.nlm.nih.gov/32666402/>

## 3. Neurosteroids and the GABA_A receptor — the central mechanism of the model

16. Maguire J, Mody I. **GABA(A)R plasticity during pregnancy: relevance to postpartum
    depression.** *Neuron* 2008. — **the key paper of the model.** Downregulation of the δ subunit in pregnancy and
    the delayed recovery after delivery; postpartum depression and abnormal maternal behaviour in the δ-KO mouse.
    <https://pubmed.ncbi.nlm.nih.gov/18667149/>
17. Maguire JL, Stell BM, Rafizadeh M, Mody I. **Ovarian cycle-linked changes in
    GABA(A) receptors mediating tonic inhibition alter seizure susceptibility and
    anxiety.** *Nat Neurosci* 2005. — tonic inhibition and the steroid-dependent plasticity of δ expression.
    <https://pubmed.ncbi.nlm.nih.gov/15895085/>
18. Maguire J, Mody I. **Neurosteroid synthesis-mediated regulation of GABA(A)
    receptors: relevance to the ovarian cycle and stress.** *J Neurosci* 2007.
    <https://pubmed.ncbi.nlm.nih.gov/17329412/>
19. Hosie AM, Wilkins ME, da Silva HMA, Smart TG. **Endogenous neurosteroids regulate
    GABAA receptors through two discrete transmembrane sites.** *Nature* 2006. —
    **the structural basis for the model separating efficacy (site 1, nM) from sedation (site 2, µM).**
    <https://pubmed.ncbi.nlm.nih.gov/17108970/>
20. Hosie AM, et al. **Neurosteroid binding sites on GABA(A) receptors.** *Pharmacol
    Ther* 2007. <https://pubmed.ncbi.nlm.nih.gov/17560657/>
21. Belelli D, Lambert JJ. **Neurosteroids: endogenous regulators of the GABA(A)
    receptor.** *Nat Rev Neurosci* 2005.
    <https://pubmed.ncbi.nlm.nih.gov/15959466/>
22. Stell BM, Brickley SG, Tang CY, Farrant M, Mody I. **Neuroactive steroids reduce
    neuronal excitability by selectively enhancing tonic inhibition mediated by delta
    subunit-containing GABAA receptors.** *PNAS* 2003. — extrasynaptic selectivity.
    <https://pubmed.ncbi.nlm.nih.gov/14623958/>
23. Mihalek RM, et al. **Attenuated sensitivity to neuroactive steroids in
    gamma-aminobutyrate type A receptor delta subunit knockout mice.** *PNAS* 1999.
    <https://pubmed.ncbi.nlm.nih.gov/10536021/>
24. Concas A, et al. **Role of brain allopregnanolone in the plasticity of
    gamma-aminobutyric acid type A receptor in rat brain during pregnancy and after
    delivery.** *PNAS* 1998. — **the experimental prototype of the model's two time constants (ligand vs receptor).**
    <https://pubmed.ncbi.nlm.nih.gov/9789080/>
25. Sarkar J, Wakefield S, MacKenzie G, Moss SJ, Maguire J. **Neurosteroidogenesis is
    required for the physiological response to stress: role of neurosteroid-sensitive
    GABAA receptors.** *J Neurosci* 2011. — the δ-GABA_A brake on PVN CRH neurons;
    the HPA-circuit coupling of the model.
    <https://pubmed.ncbi.nlm.nih.gov/22171026/>
26. Ferando I, Mody I. **Interneuronal GABAA receptors inside and outside of
    synapses.** *Curr Opin Neurobiol* 2014.
    <https://pubmed.ncbi.nlm.nih.gov/24650505/>
27. Paul SM, Purdy RH. **Neuroactive steroids.** *FASEB J* 1992.
    <https://pubmed.ncbi.nlm.nih.gov/1347506/>
28. Luscher B, Shen Q, Sahir N. **The GABAergic deficit hypothesis of major depressive
    disorder.** *Mol Psychiatry* 2011.
    <https://pubmed.ncbi.nlm.nih.gov/21079608/>
29. A review that includes Windle RJ and others: **[GABAergic approach of postpartum depression: A
    translational review of literature].** *Encephale* 2020. — a translational synthesis that includes KCC2 and
    chloride homeostasis. <https://pubmed.ncbi.nlm.nih.gov/31767256/>
30. **The human gamma-aminobutyric acid A receptor delta (GABRD) gene: molecular
    characterisation and tissue-specific expression.** *Gene* 2002. — the δ subunit gene.
    <https://pubmed.ncbi.nlm.nih.gov/12119096/>
31. Bäckström T, et al. **Positive GABA(A) receptor modulating steroids and their
    antagonists: implications for clinical treatments.** *J Neuroendocrinol* 2022. —
    the functional antagonism of isoallopregnanolone (sepranolone); the ISOALLO node of the model.
    <https://pubmed.ncbi.nlm.nih.gov/34337790/>
32. Bixo M, et al. **Treatment of premenstrual dysphoric disorder with the GABA(A)
    receptor modulating steroid antagonist Sepranolone (UC1010).**
    *Psychoneuroendocrinology* 2017.
    <https://pubmed.ncbi.nlm.nih.gov/28319848/>

## 4. Steroid metabolism and synthesis pathways — the synthesis and clearance parameters

33. **Serum concentrations of progesterone and 5 alpha-pregnane-3,20-dione during
    labor and early post partum.** *Acta Obstet Gynecol Scand* 1990. — the basis for the elimination half-lives
    of P4 and 5α-DHP in the model.
    <https://pubmed.ncbi.nlm.nih.gov/2386015/>
34. **Allopregnanolone serum levels and expression of 5 alpha-reductase and 3
    alpha-hydroxysteroid dehydrogenase isoforms.** *Epilepsy Res* 2003.
    <https://pubmed.ncbi.nlm.nih.gov/12742591/>
35. **Hypothalamic 5 alpha-reductase and 3 alpha-oxidoreductase activity.**
    *J Steroid Biochem Mol Biol* 2002. — local synthesis in the brain (the model's BRSYN term).
    <https://pubmed.ncbi.nlm.nih.gov/11867268/>
36. Timby E, et al. **Pharmacokinetic and behavioral effects of allopregnanolone in
    healthy women.** *Psychopharmacology (Berl)* 2006. — the PK of IV ALLO in humans and its
    sedative effect. <https://pubmed.ncbi.nlm.nih.gov/16177884/>
37. **Allopregnanolone assays.** *J Clin Endocrinol Metab* 2001. — RIA vs LC-MS/MS
    discrepancy; the basis for the model's assumption A4 (why the ratio rather than the absolute value was matched).
    <https://pubmed.ncbi.nlm.nih.gov/11238543/>
38. **Peripartum allopregnanolone blood concentrations and depressive symptoms: a
    systematic review and individual participant data meta-analysis.** *Mol
    Psychiatry* 2025. — **the key point:** the relationship between the absolute ALLO concentration and symptoms is not consistent →
    supporting the model's premise that the problem is not "low ALLO" but the rate of change and the receptor response.
    <https://pubmed.ncbi.nlm.nih.gov/39511449/>
39. **Allopregnanolone in the peripartum: correlates, concentrations, and challenges —
    a systematic review.** *Psychoneuroendocrinology* 2024.
    <https://pubmed.ncbi.nlm.nih.gov/38759520/>
40. **Profiling neuroactive steroids in pregnancy: a non-derivatised LC-MS/MS
    method.** *J Chromatogr B* 2025.
    <https://pubmed.ncbi.nlm.nih.gov/40054418/>

## 5. Human data — hormone withdrawal and vulnerability

41. Bloch M, et al. **Effects of gonadal steroids in women with a history of postpartum
    depression.** *Am J Psychiatry* 2000. — **the experimental basis for the model's concept of vulnerability (V):**
    the same steroid load and withdrawal provokes symptoms only in women with a history.
    <https://pubmed.ncbi.nlm.nih.gov/10831472/>
42. Bloch M, Daly RC, Rubinow DR. **Endocrine factors in the etiology of postpartum
    depression.** *Compr Psychiatry* 2003.
    <https://pubmed.ncbi.nlm.nih.gov/12764712/>
43. Schiller CE, Meltzer-Brody S, Rubinow DR. **The role of reproductive hormones in
    postpartum depression.** *CNS Spectr* 2015.
    <https://pubmed.ncbi.nlm.nih.gov/25263255/>
44. Nappi RE, et al. **Serum allopregnanolone in women with postpartum "blues".**
    *Obstet Gynecol* 2001. <https://pubmed.ncbi.nlm.nih.gov/11152912/>
45. Hellgren C, et al. **Low serum allopregnanolone is associated with symptoms of
    depression in late pregnancy.** *Neuropsychobiology* 2014.
    <https://pubmed.ncbi.nlm.nih.gov/24776841/>
46. Osborne LM, et al. **Lower allopregnanolone during pregnancy predicts postpartum
    depression: an exploratory study.** *Psychoneuroendocrinology* 2017.
    <https://pubmed.ncbi.nlm.nih.gov/28278440/>
47. Epperson CN, et al. **Preliminary evidence of reduced occipital GABA
    concentrations in puerperal women: a 1H-MRS study.** *Psychopharmacology (Berl)*
    2006. <https://pubmed.ncbi.nlm.nih.gov/16724188/>
48. Deligiannidis KM, et al. **Resting-state functional connectivity, cortical GABA,
    and neuroactive steroids in peripartum and peripartum depressed women.**
    *Neuropsychopharmacology* 2019.
    <https://pubmed.ncbi.nlm.nih.gov/30327498/>
49. Deligiannidis KM, et al. **GABAergic neuroactive steroids and resting-state
    functional connectivity in postpartum depression: a preliminary study.**
    *J Psychiatr Res* 2013.
    <https://pubmed.ncbi.nlm.nih.gov/23499388/>

## 6. The HPA axis (the second slow variable)

50. Magiakou MA, et al. **Hypothalamic corticotropin-releasing hormone suppression
    during the postpartum period: implications for the increase in psychiatric
    manifestations at this time.** *J Clin Endocrinol Metab* 1996. — **the basis for the model's HCRH
    recovery time constant (a window of ≈ 12 weeks).**
    <https://pubmed.ncbi.nlm.nih.gov/8626857/>
51. Yim IS, et al. **Risk of postpartum depressive symptoms with elevated
    corticotropin-releasing hormone in human pregnancy.** *Arch Gen Psychiatry* 2009.
    <https://pubmed.ncbi.nlm.nih.gov/19188538/>
52. Glynn LM, Davis EP, Sandman CA. **New insights into the role of perinatal HPA-axis
    dysregulation in postpartum depression.** *Neuropeptides* 2013.
    <https://pubmed.ncbi.nlm.nih.gov/24210135/>

## 7. Monoamines · inflammation · kynurenine

53. Sacher J, et al. **Elevated brain monoamine oxidase A binding in the early
    postpartum period.** *Arch Gen Psychiatry* 2010. — **the direct basis for the model's MAO-A +21 %
    term** (E2 withdrawal → a rise in MAO-A).
    <https://pubmed.ncbi.nlm.nih.gov/20439828/>
54. Osborne LM, Monk C. **Perinatal depression — the fourth inflammatory morbidity of
    pregnancy? Theory and literature review.** *Psychoneuroendocrinology* 2013.
    <https://pubmed.ncbi.nlm.nih.gov/23608136/>
55. Bränn E, et al. **Inflammatory markers in women with postpartum depressive
    symptoms.** *J Neurosci Res* 2020.
    <https://pubmed.ncbi.nlm.nih.gov/30252150/>
56. Achtyes E, et al. **Inflammation and kynurenine pathway dysregulation in
    post-partum women with severe and suicidal depression.** *Brain Behav Immun*
    2020. — the model's KYN/TRP and QUIN/KYNA terms.
    <https://pubmed.ncbi.nlm.nih.gov/31698012/>
57. **Associations between estrogen and progesterone, the kynurenine pathway, and
    inflammation in the post-partum.** *J Affect Disord* 2021.
    <https://pubmed.ncbi.nlm.nih.gov/33278766/>
58. Corwin EJ, et al. **Symptoms of postpartum depression associated with elevated
    levels of interleukin-1 beta during the first month postpartum.** *Biol Res Nurs*
    2008. <https://pubmed.ncbi.nlm.nih.gov/18829596/>
59. Corwin EJ, Pajer K. **The psychoneuroimmunology of postpartum depression.**
    *J Womens Health* 2008.
    <https://pubmed.ncbi.nlm.nih.gov/18710365/>

## 8. Genetics and epigenetics (vulnerability, V)

60. Guintivano J, Arad M, Gould TD, Payne JL, Kaminsky ZA. **Antenatal prediction of
    postpartum depression with blood DNA methylation biomarkers.** *Mol Psychiatry*
    2014. — HP1BP3 / TTC9B.
    <https://pubmed.ncbi.nlm.nih.gov/23689534/>
61. **Meta-Analyses of Genome-Wide Association Studies for Postpartum Depression.**
    *Am J Psychiatry* 2023. — the first large-scale GWAS meta-analysis of PPD.
    <https://pubmed.ncbi.nlm.nih.gov/37849304/>
62. **The First Large GWAS Meta-Analysis for Postpartum Depression** (editorial).
    *Am J Psychiatry* 2023.
    <https://pubmed.ncbi.nlm.nih.gov/38037399/>
63. Viktorin A, et al. **Heritability of perinatal depression and genetic overlap with
    nonperinatal depression.** *Am J Psychiatry* 2016.
    <https://pubmed.ncbi.nlm.nih.gov/26337037/>

## 9. Sleep — the vicious circle of the model (the closed loop)

64. Lawson A, et al. **The relationship between sleep and postpartum mental disorders:
    a systematic review.** *J Affect Disord* 2015.
    <https://pubmed.ncbi.nlm.nih.gov/25702602/>
65. **Changes in sleep quality, but not hormones predict time to postpartum depression
    recurrence.** *J Affect Disord* 2011. — **the direct basis for the model placing sleep, rather than
    hormones, as the strongest proximal driver.**
    <https://pubmed.ncbi.nlm.nih.gov/20708275/>
66. Okun ML, et al. **Sleep complaints in late pregnancy and the recurrence of
    postpartum depression.** *Behav Sleep Med* 2009.
    <https://pubmed.ncbi.nlm.nih.gov/19330583/>
67. Goyal D, Gay C, Lee K. **Fragmented maternal sleep is more strongly correlated with
    depressive symptoms than infant temperament at three months postpartum.** *Arch
    Womens Ment Health* 2009. — fragmentation matters more than total sleep time.
    <https://pubmed.ncbi.nlm.nih.gov/19396527/>
68. **Maternal depressive symptoms, dysfunctional cognitions, and infant night waking:
    the role of maternal nighttime behavior.** *Child Dev* 2012. — a bidirectional loop.
    <https://pubmed.ncbi.nlm.nih.gov/22506917/>
69. **Preventing recurrence of postpartum depression by regulating sleep.** *Expert Rev
    Neurother* 2023.
    <https://pubmed.ncbi.nlm.nih.gov/37462620/>

## 10. Oxytocin · mother-infant bonding · infant outcomes

70. Skrundz M, et al. **Plasma oxytocin concentration during pregnancy is associated
    with development of postpartum depression.** *Neuropsychopharmacology* 2011.
    <https://pubmed.ncbi.nlm.nih.gov/21562482/>
71. Moehler E, et al. **Maternal depressive symptoms in the postnatal period are
    associated with long-term impairment of mother-infant bonding.** *Arch Womens Ment
    Health* 2006. <https://pubmed.ncbi.nlm.nih.gov/16937313/>
72. **Associations between maternal psychological distress and mother-infant bonding: a
    systematic review and meta-analysis.** *Arch Womens Ment Health* 2023.
    <https://pubmed.ncbi.nlm.nih.gov/37316760/>
73. Netsi E, et al. **Association of persistent and severe postnatal depression with
    child outcomes.** *JAMA Psychiatry* 2018. — the outcomes in the children of persistent, severe PPD.
    <https://pubmed.ncbi.nlm.nih.gov/29387878/>

## 11. Neurosteroid therapeutic trials — the model's calibration targets

74. Kanes S, et al. **Brexanolone (SAGE-547 injection) in post-partum depression: a
    randomised controlled trial.** *Lancet* 2017. — phase 2, n = 21, severe
    (HAM-D ≥ 26): **change in HAM-D at 60 hours −21.0 (brexanolone) vs −8.8 (placebo).**
    <https://pubmed.ncbi.nlm.nih.gov/28619476/>
75. Meltzer-Brody S, et al. **Brexanolone injection in post-partum depression: two
    multicentre, double-blind, randomised, placebo-controlled, phase 3 trials.**
    *Lancet* 2018. — **study 1 (HAM-D ≥ 26): at 60 hours −19.5 (BRX60) · −17.7 (BRX90) ·
    −14.0 (placebo); study 2 (HAM-D 20-25): −14.6 vs −12.1.** The observation that the dose response from
    60 to 90 µg/kg/h is flat corresponds to the Hill saturation the model predicts.
    <https://pubmed.ncbi.nlm.nih.gov/30177236/>
75b. Kanes SJ, et al. **Open-label, proof-of-concept study of brexanolone in the
    treatment of severe postpartum depression.** *Hum Psychopharmacol* 2017.
    <https://pubmed.ncbi.nlm.nih.gov/28370307/>
76. Deligiannidis KM, et al. **Effect of zuranolone vs placebo in postpartum
    depression: a randomized clinical trial (ROBIN).** *JAMA Psychiatry* 2021. —
    zuranolone 30 mg × 14 days: **day 15 HAM-D −17.8 vs −13.6.**
    <https://pubmed.ncbi.nlm.nih.gov/34190962/>
77. Deligiannidis KM, et al. **Zuranolone for the Treatment of Postpartum Depression
    (SKYLARK).** *Am J Psychiatry* 2023. — zuranolone 50 mg × 14 days: **day 15 HAM-D
    −15.6 vs −11.6 (a difference of −4.0); significant at days 3, 28, and 45 as well.**
    <https://pubmed.ncbi.nlm.nih.gov/37491938/>
78. Althaus AL, et al. **Preclinical characterization of zuranolone (SAGE-217), a
    selective neuroactive steroid GABA(A) receptor positive allosteric modulator.**
    *Neuropharmacology* 2020. — simultaneous synaptic and extrasynaptic action (the model's W_PHASIC term).
    <https://pubmed.ncbi.nlm.nih.gov/32976892/>
79. Gunduz-Bruce H, et al. **Trial of SAGE-217 in patients with major depressive
    disorder.** *N Engl J Med* 2019.
    <https://pubmed.ncbi.nlm.nih.gov/31483961/>
80. **Efficacy and safety of zuranolone co-initiated with an antidepressant in adults
    with major depressive disorder (CORAL).** *Neuropsychopharmacology* 2024.
    <https://pubmed.ncbi.nlm.nih.gov/37875578/>
81. **Zuranolone for treatment of major depressive disorder: a systematic review and
    meta-analysis.** *Front Neurosci* 2024.
    <https://pubmed.ncbi.nlm.nih.gov/38726035/>
82. **Brexanolone, zuranolone and related neurosteroid GABA(A) receptor positive
    allosteric modulators for postnatal depression.** *Cochrane Database Syst Rev*
    2025. — the most recent systematic synthesis.
    <https://pubmed.ncbi.nlm.nih.gov/40562419/>
83. **Preclinical and clinical pharmacology of brexanolone (allopregnanolone) for
    postpartum depression: a landmark journey.** *Psychopharmacology (Berl)* 2023.
    <https://pubmed.ncbi.nlm.nih.gov/37566239/>
84. Cornett EM, et al. **Brexanolone to treat postpartum depression in adult women.**
    *Psychopharmacol Bull* 2021.
    <https://pubmed.ncbi.nlm.nih.gov/34092826/>
85. Faden J, Citrome L. **Brexanolone for postpartum depression: clinical evidence and
    practical considerations.** *Pharmacotherapy* 2019. — the practicalities of REMS and of monitoring for sedation.
    <https://pubmed.ncbi.nlm.nih.gov/31514247/>
86. **Pharmacokinetics, safety, and tolerability of single and multiple doses of
    zuranolone in Japanese and White healthy subjects.** *Neuropsychopharmacol Rep*
    2023. — the basis for the model's assumption that zuranolone t½ ≈ 20 h.
    <https://pubmed.ncbi.nlm.nih.gov/37366077/>
87. **The magnitude and sustainability of treatment benefit of zuranolone on function
    and well-being (SDS).** *J Affect Disord* 2024.
    <https://pubmed.ncbi.nlm.nih.gov/38325605/>
88. **Real-world safety profile of zuranolone for postpartum depression: a FAERS
    analysis.** *J Affect Disord* 2026.
    <https://pubmed.ncbi.nlm.nih.gov/40935249/>

## 12. Drug transfer in lactation — the infant-exposure arm

89. **Allopregnanolone Concentrations in Breast Milk and Plasma from Healthy Volunteers
    Receiving Brexanolone Injection, With Population Pharmacokinetic Modeling of
    Potential Relative Infant Dose.** *Clin Pharmacokinet* 2022. — **measured data directly comparable with
    the RID calculation of the model.**
    <https://pubmed.ncbi.nlm.nih.gov/35869362/>
90. **Using Brexanolone for Postpartum Depression Must Account for Lactation.**
    *Matern Child Health J* 2021.
    <https://pubmed.ncbi.nlm.nih.gov/34019187/>
91. Weissman AM, et al. **Pooled analysis of antidepressant levels in lactating
    mothers, breast milk, and nursing infants.** *Am J Psychiatry* 2004. — sertraline has
    the lowest RID among the SSRIs.
    <https://pubmed.ncbi.nlm.nih.gov/15169695/>
92. Berle JØ, et al. **Breastfeeding during maternal antidepressant treatment with
    serotonin reuptake inhibitors: infant exposure, clinical symptoms, and cytochrome
    P450 genotypes.** *J Clin Psychiatry* 2004.
    <https://pubmed.ncbi.nlm.nih.gov/15367050/>

## 13. Existing and non-pharmacological treatments (comparator arms)

93. Molyneaux E, et al. **Antidepressant treatment for postnatal depression.**
    *Cochrane Database Syst Rev* 2014.
    <https://pubmed.ncbi.nlm.nih.gov/25211400/>
94. Hantsoo L, et al. **A randomized, placebo-controlled, double-blind trial of
    sertraline for postpartum depression.** *Psychopharmacology (Berl)* 2014.
    <https://pubmed.ncbi.nlm.nih.gov/24173623/>
95. Wisner KL, et al. **Prevention of postpartum depression: a pilot randomized
    clinical trial.** *Am J Psychiatry* 2004.
    <https://pubmed.ncbi.nlm.nih.gov/15229064/>
96. Gregoire AJ, et al. **Transdermal oestrogen for treatment of severe postnatal
    depression.** *Lancet* 1996. — the E2 patch scenario of the model.
    <https://pubmed.ncbi.nlm.nih.gov/8598756/>
97. **Perioperative Adjunctive Esketamine for Postpartum Depression Among Women
    Undergoing Elective Cesarean Delivery: A Randomized Clinical Trial.** *JAMA Netw
    Open* 2024. <https://pubmed.ncbi.nlm.nih.gov/38446480/>
98. **Intraoperative Esketamine and Postpartum Depression Among Women With Cesarean
    Delivery: A Randomized Clinical Trial.** *JAMA Netw Open* 2025.
    <https://pubmed.ncbi.nlm.nih.gov/39946130/>
99. **Prophylactic esketamine for postpartum depression after cesarean section: a
    systematic review and meta-analysis.** *J Affect Disord* 2025.
    <https://pubmed.ncbi.nlm.nih.gov/40505985/>
100. O'Connor E, et al. **Interventions to Prevent Perinatal Depression: Evidence
     Report and Systematic Review for the US Preventive Services Task Force.** *JAMA*
     2019. <https://pubmed.ncbi.nlm.nih.gov/30747970/>
101. Dennis CL, Dowswell T. **Psychosocial and psychological interventions for
     preventing postpartum depression.** *Cochrane Database Syst Rev* 2013.
     <https://pubmed.ncbi.nlm.nih.gov/23450532/>
102. Zlotnick C, et al. **Postpartum depression in women receiving public assistance:
     pilot study of an interpersonal-therapy-oriented group intervention (ROSE).**
     *Am J Psychiatry* 2001.
     <https://pubmed.ncbi.nlm.nih.gov/11282702/>
103. **Effect of Bright Light Therapy on Perinatal Depression: A Systematic Review and
     Meta-Analysis.** *Can J Psychiatry* 2024.
     <https://pubmed.ncbi.nlm.nih.gov/38863243/>
104. **The effectiveness of exercise-based interventions for preventing or treating
     postpartum depression: a systematic review and meta-analysis.** *Arch Womens Ment
     Health* 2019. <https://pubmed.ncbi.nlm.nih.gov/29882074/>
105. **Effectiveness of aerobic exercise in the prevention and treatment of postpartum
     depression: meta-analysis and network meta-analysis.** *PLoS One* 2023.
     <https://pubmed.ncbi.nlm.nih.gov/38019729/>
106. Muller AF, et al. **Postpartum thyroiditis.** *Best Pract Res Clin Endocrinol
     Metab* 2004. — differential diagnosis and comorbidity (the model's THYRPP node).
     <https://pubmed.ncbi.nlm.nih.gov/15157842/>

## 14. Tools and methodology

107. gPKPDviz — an mrgsolve-based Shiny tool for PK/PD simulation.
     <https://pmc.ncbi.nlm.nih.gov/articles/PMC10941578/>
108. mrgsolve — R package for simulation from ODE-based models.
     <https://mrgsolve.org/>
109. **New directions in neurosteroid therapeutics in neuropsychiatry.** *Neurosci
     Biobehav Rev* 2025. — a prospect for the next generation of candidate compounds.
     <https://pubmed.ncbi.nlm.nih.gov/40127877/>

---

## What the model does **not** get from the literature

Stated here for the sake of honesty.

| Item | Status |
|------|------|
| `KR` (δ receptor recovery t½ ≈ 7 days) | **No** human data. A free parameter inferred from the rodent (#24). The most important parameter of the model, and first in the sensitivity analysis. |
| `THR0` (the E/I symptom threshold) | An unobservable latent parameter. Calibrated to the placebo arm trajectory. |
| `ZUR_EQ` (zuranolone → ALLO equivalence) | The **single** calibration parameter, fitted to SKYLARK day 15. |
| `K_CARE`, `K_NSP` (non-specific effects) | Phenomenological terms fitted to the placebo arm rather than to a mechanism. Held fixed in the active arms. |
| `W_SELF` (symptom self-reinforcement) | A shorthand bundling rumination, the DMN, and glutamate into one. Each has its own evidence but the coefficient is free. |
| Hypothalamic CRH recovery over 12 weeks | Rests on a single small study, #50. |
| F_oral = 0.05 for infant exposure | No neonatal measurement. A conservative assumption, comparable only against the measured RID of #89. |
