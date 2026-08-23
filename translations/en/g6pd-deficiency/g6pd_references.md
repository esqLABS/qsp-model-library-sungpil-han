# G6PD deficiency QSP model — references
# Glucose-6-Phosphate Dehydrogenase Deficiency — Reference List

The sources of every structural claim and every parameter central value in this model. Section 1 matters most —
the whole model rests on one line: "enzyme activity is not a number but **a function of cell age**, and a
variant is the **slope** of that function" — and the paper that actually measured that slope in
humans is in section 1.

Every structural claim and central parameter value in
`g6pd_mrgsolve_model.R` traces to something here. Section 1 carries the
model's entire thesis: the papers that measured, *in vivo*, the rate at
which G6PD decays inside a red cell that can no longer make any.

---

## 1. The model's central claim — the enzyme is a function of cell age (THE THESIS)

**The first two papers of this section are the model's `E0` and `TAU` parameters themselves.**

1. **Piomelli S, Corash LM, Davenport DD, Miraglia J, Amorosi EL.**
   In vivo lability of glucose-6-phosphate dehydrogenase in GdA− and
   Gd Mediterranean deficiency. *J Clin Invest.* 1968;47(4):940–948.
   <https://pubmed.ncbi.nlm.nih.gov/5641628/>
   — This model's `TAU`. Red cells were separated by age on a density gradient and their enzyme
   activity measured, showing directly a half-life of about 62 days for normal type B against about 13 days for A−. A mature red
   cell cannot make new protein, so the enzyme can only decline, and activity therefore becomes
   an exponential function of age. **The model's single line `E(a) = E0·exp(−ln2·a/TAU)`
   comes entirely from this paper.**

2. **Dern RJ, Beutler E, Alving AS.** The hemolytic effect of primaquine.
   II. The natural course of the hemolytic anemia and the mechanism of its
   self-limited character. *J Lab Clin Med.* 1954;44(2):171–176.
   <https://pubmed.ncbi.nlm.nih.gov/13184228/>
   — The original of scenario 1. Primaquine 30 mg/day was continued **without being stopped**, and yet
   the haemoglobin bottomed out at 7–8 days and **recovered while the drug was still being taken**. As early as
   1954 the authors named the reason exactly: the surviving cells are young, and the newly released
   reticulocytes are full of enzyme. The model reproduces this curve without a single parameter
   that knows about the phenomenon.

3. **Beutler E.** Glucose-6-phosphate dehydrogenase deficiency: a historical
   perspective. *Blood.* 2008;111(1):16–24.
   <https://pubmed.ncbi.nlm.nih.gov/18156501/>

4. **Morelli A, Benatti U, Gaetani GF, De Flora A.** Biochemical mechanisms
   of glucose-6-phosphate dehydrogenase deficiency. *Proc Natl Acad Sci USA.*
   1978;75(4):1979–1983.
   <https://pubmed.ncbi.nlm.nih.gov/273925/>
   — The paper that separated the **instability** of a variant enzyme (a fast TAU) from its **low initial activity**
   (a low E0) and showed them to be distinct biochemical defects. Why the model keeps the two
   parameters separate.

5. **Alving AS, Johnson CF, Tarlov AR, Brewer GJ, Kellermeyer RW, Carson PE.**
   Mitigation of the haemolytic effect of primaquine and enhancement of its
   action against exoerythrocytic forms of the Chesson strain of *Plasmodium
   vivax* by intermittent regimens of drug administration.
   *Bull World Health Organ.* 1960;22:621–631.
   <https://pubmed.ncbi.nlm.nih.gov/13793053/>
   — The basis for scenario 5 (once-weekly dosing). That weekly dosing is not simply a dose reduction but
   **a different mechanism** — each pulse shaves off only the thinnest, oldest layer, and the 7-day interval gives
   the reticulocytes time to fill that place.

---

## 2. Reviews · the whole clinical picture

6. **Luzzatto L, Ally M, Notaro R.** Glucose-6-phosphate dehydrogenase
   deficiency. *Blood.* 2020;136(11):1225–1240.
   <https://pubmed.ncbi.nlm.nih.gov/32702756/>
   — The current standard review.

7. **Luzzatto L, Arese P.** Favism and glucose-6-phosphate dehydrogenase
   deficiency. *N Engl J Med.* 2018;378(1):60–71.
   <https://pubmed.ncbi.nlm.nih.gov/29298156/>
   — The clinical description and time course of scenario 6 (favism) (a few hours to 48 hours).

8. **Cappellini MD, Fiorelli G.** Glucose-6-phosphate dehydrogenase
   deficiency. *Lancet.* 2008;371(9606):64–74.
   <https://pubmed.ncbi.nlm.nih.gov/18177777/>

9. **Beutler E.** G6PD deficiency. *Blood.* 1994;84(11):3613–3636.
   <https://pubmed.ncbi.nlm.nih.gov/7949118/>

10. **Bancone G, Chu CS.** G6PD variants and haemolytic sensitivity to
    primaquine and other drugs. *Front Pharmacol.* 2021;12:638885.
    <https://pubmed.ncbi.nlm.nih.gov/33790795/>
    — Addresses head-on the point that haemolytic susceptibility by variant is not predicted from
    "% activity" alone. The same claim as the model's A− versus Mediterranean contrast (scenario 1 vs 2).

---

## 3. Genetics · classification of variants · epidemiology

11. **WHO Working Group.** Glucose-6-phosphate dehydrogenase deficiency.
    *Bull World Health Organ.* 1989;67(6):601–611.
    <https://pubmed.ncbi.nlm.nih.gov/2633878/>
    — The classic Class I–V classification.

12. **Luzzatto L, Bancone G, Dugué P-A, et al.** New WHO classification of
    genetic variants causing G6PD deficiency.
    *Bull World Health Organ.* 2024;102(8):615–617.
    <https://pubmed.ncbi.nlm.nih.gov/39070596/>
    — The 2022 revised classification (on a % activity basis). The framework the model's `VARIANTS` list follows.

13. **Howes RE, Piel FB, Patil AP, et al.** G6PD deficiency prevalence and
    estimated affected population size: a geostatistical model-based map and
    population estimates. *PLoS Med.* 2012;9(11):e1001339.
    <https://pubmed.ncbi.nlm.nih.gov/23152723/>
    — About 400 million people worldwide.

14. **Nkhoma ET, Poole C, Vannappagari V, Hall SA, Beutler E.** The global
    prevalence of glucose-6-phosphate dehydrogenase deficiency: a systematic
    review and meta-analysis. *Blood Cells Mol Dis.* 2009;42(3):267–278.
    <https://pubmed.ncbi.nlm.nih.gov/19233695/>

15. **Minucci A, Moradkhani K, Hwang MJ, Zuppi C, Giardina B, Capoluongo E.**
    Glucose-6-phosphate dehydrogenase (G6PD) mutations database: review of
    the "old" and update of the new mutations.
    *Blood Cells Mol Dis.* 2012;48(3):154–165.
    <https://pubmed.ncbi.nlm.nih.gov/22293322/>

16. **Beutler E, Vulliamy TJ.** Hematologically important mutations:
    glucose-6-phosphate dehydrogenase.
    *Blood Cells Mol Dis.* 2002;28(2):93–103.
    <https://pubmed.ncbi.nlm.nih.gov/12064901/>

17. **Ruwende C, Khoo SC, Snow RW, et al.** Natural selection of hemi- and
    heterozygotes for G6PD deficiency in Africa by resistance to severe
    malaria. *Nature.* 1995;376(6537):246–249.
    <https://pubmed.ncbi.nlm.nih.gov/7617034/>
    — Why this allele is so common.

---

## 4. Heterozygous females · lyonisation

**Why the model's `mix_mosaic()` exists.** A heterozygous female is not one phenotype
but **two red cell populations**, and a whole-blood assay averages the two.

18. **Chu CS, Bancone G, Nosten F, White NJ, Luzzatto L.** Primaquine-induced
    haemolysis in females heterozygous for G6PD deficiency.
    *Malar J.* 2018;17(1):101.
    <https://pubmed.ncbi.nlm.nih.gov/29499733/>

19. **Chu CS, Bancone G, Moore KA, et al.** Haemolysis in G6PD heterozygous
    females treated with primaquine for *Plasmodium vivax* malaria: a nested
    cohort in a trial of radical curative regimens.
    *PLoS Med.* 2017;14(2):e1002224.
    <https://pubmed.ncbi.nlm.nih.gov/28170391/>

20. **Rueangweerayut R, Bancone G, Harrell EJ, et al.** Hemolytic potential of
    tafenoquine in female volunteers heterozygous for glucose-6-phosphate
    dehydrogenase (G6PD) deficiency (G6PD Mahidol variant) versus G6PD-normal
    volunteers. *Am J Trop Med Hyg.* 2017;97(3):702–711.
    <https://pubmed.ncbi.nlm.nih.gov/28749773/>
    — The clinical counterpart of scenario 11. A fall in haemoglobin was observed after a single 300 mg dose
    of tafenoquine in heterozygous females with 40–60% activity.

21. **Bancone G, Kalnoky M, Chu CS, et al.** The G6PD flow-cytometric assay is
    a reliable tool for diagnosis of G6PD deficiency in women and anaemic
    subjects. *Sci Rep.* 2017;7(1):9822.
    <https://pubmed.ncbi.nlm.nih.gov/28852037/>
    — The only assay that sees the **distribution**. Direct evidence for why a quantitative assay that sees only
    the mean fails in women.

---

## 5. The biochemistry of oxidative damage — GSH · haemichrome · band 3

22. **Arese P, Gallo V, Pantaleo A, Turrini F.** Life and death of
    glucose-6-phosphate dehydrogenase (G6PD) deficient erythrocytes — role of
    redox stress and band 3 modifications.
    *Transfus Med Hemother.* 2012;39(5):328–334.
    <https://pubmed.ncbi.nlm.nih.gov/23801924/>
    — The skeleton of model clusters 6–8: GSH depletion → haemichrome → Heinz body →
    band 3 clustering → autologous IgG/complement → macrophage clearance.

23. **Turrini F, Arese P, Yuan J, Low PS.** Clustering of integral membrane
    proteins of the human erythrocyte membrane stimulates autologous IgG
    binding, complement deposition, and phagocytosis.
    *J Biol Chem.* 1991;266(35):23611–23617.
    <https://pubmed.ncbi.nlm.nih.gov/1748639/>

24. **Low FM, Hampton MB, Winterbourn CC.** Peroxiredoxin 2 and peroxide
    metabolism in the erythrocyte.
    *Antioxid Redox Signal.* 2008;10(9):1621–1630.
    <https://pubmed.ncbi.nlm.nih.gov/18479200/>
    — The order of magnitude of the model's `VMAXOX` (the peroxide-handling capacity of a normal young cell).

25. **Gaetani GF, Ferraris AM, Rolfo M, Mangerini R, Arena S, Kirkman HN.**
    Predominant role of catalase in the disposal of hydrogen peroxide within
    human erythrocytes. *Blood.* 1996;87(4):1595–1599.
    <https://pubmed.ncbi.nlm.nih.gov/8608252/>

26. **Rifkind JM, Nagababu E.** Hemoglobin redox reactions and red blood cell
    aging. *Antioxid Redox Signal.* 2013;18(17):2274–2283.
    <https://pubmed.ncbi.nlm.nih.gov/23025272/>
    — The model's `OXBASE` (the baseline autoxidation load present even with no drug).

27. **Lang F, Lang KS, Lang PA, Huber SM, Wieder T.** Mechanisms and
    significance of eryptosis. *Antioxid Redox Signal.* 2006;8(7–8):1183–1192.
    <https://pubmed.ncbi.nlm.nih.gov/16910766/>

28. **Fibach E, Rachmilewitz E.** The role of oxidative stress in hemolytic
    anemia. *Curr Mol Med.* 2008;8(7):609–619.
    <https://pubmed.ncbi.nlm.nih.gov/18991647/>

29. **Tang HY, Ho HY, Wu PR, et al.** Inability to maintain GSH pool in G6PD-
    deficient red cells causes futile AMPK activation and irreversible
    metabolic disturbance. *Antioxid Redox Signal.* 2015;22(9):744–759.
    <https://pubmed.ncbi.nlm.nih.gov/25580850/>

---

## 6. Primaquine · tafenoquine (8-aminoquinolines: the half-life IS the toxicology)

30. **Watson J, Taylor WRJ, Menard D, Kheng S, White NJ.** Modelling
    primaquine-induced haemolysis in G6PD deficiency.
    *eLife.* 2017;6:e23061.
    <https://pubmed.ncbi.nlm.nih.gov/28553934/>
    — **The prior quantitative model.** It too treats red cell age structure explicitly, and explains why weekly
    dosing is safe in the same way. This model amounts to adding drug
    PK (multiple drugs), the dual methaemoglobin/NADPH consumption, the bilirubin-kernicterus axis and the
    intravascular haemolysis-AKI axis on top of it.

31. **Baird JK.** 8-aminoquinoline therapy for latent malaria.
    *Clin Microbiol Rev.* 2019;32(4):e00011-19.
    <https://pubmed.ncbi.nlm.nih.gov/31366609/>

32. **Ashley EA, Recht J, White NJ.** Primaquine: the risks and the benefits.
    *Malar J.* 2014;13:418.
    <https://pubmed.ncbi.nlm.nih.gov/25363455/>

33. **Llanos-Cuentas A, Lacerda MVG, Hien TT, et al.** Tafenoquine versus
    primaquine to prevent relapse of *Plasmodium vivax* malaria.
    *N Engl J Med.* 2019;380(3):229–241.
    <https://pubmed.ncbi.nlm.nih.gov/30650326/>

34. **Lacerda MVG, Llanos-Cuentas A, Krudsood S, et al.** Single-dose
    tafenoquine to prevent relapse of *Plasmodium vivax* malaria.
    *N Engl J Med.* 2019;380(3):215–228.
    <https://pubmed.ncbi.nlm.nih.gov/30650322/>
    — The clinical context of scenario 4a. Tafenoquine 300 mg as a single dose, half-life about 15 days.

35. **Kheng S, Muth S, Taylor WRJ, et al.** Tolerability and safety of weekly
    primaquine against relapse of *Plasmodium vivax* in Cambodians with
    glucose-6-phosphate dehydrogenase deficiency.
    *BMC Med.* 2015;13:203.
    <https://pubmed.ncbi.nlm.nih.gov/26303162/>
    — The prospective clinical evidence for scenario 5.

36. **Pybus BS, Marcsisin SR, Jin X, et al.** The metabolism of primaquine to
    its active metabolite is dependent on CYP 2D6.
    *Malar J.* 2013;12:212.
    <https://pubmed.ncbi.nlm.nih.gov/23782898/>
    — The model's `AS2D6`. Because the same enzyme produces **both the therapeutic effect and the haemolysis**,
    a CYP2D6 poor metaboliser haemolyses less and also fails radical cure.

37. **Bennett JW, Pybus BS, Yadava A, et al.** Primaquine failure and
    cytochrome P-450 2D6 in *Plasmodium vivax* malaria.
    *N Engl J Med.* 2013;369(14):1381–1382.
    <https://pubmed.ncbi.nlm.nih.gov/24088113/>

38. **Mihaly GW, Ward SA, Edwards G, Orme ML, Breckenridge AM.**
    Pharmacokinetics of primaquine in man: identification of the carboxylic
    acid derivative as a major plasma metabolite.
    *Br J Clin Pharmacol.* 1984;17(4):441–446.
    <https://pubmed.ncbi.nlm.nih.gov/6721990/>
    — The model's primaquine PK parameters (`VPQ`, `CLPQ`, half-life about 7 hours).

39. **Charles BG, Miller AK, Nasveld PE, Reid MG, Harris IE, Edstein MD.**
    Population pharmacokinetics of tafenoquine during malaria prophylaxis in
    healthy subjects. *Antimicrob Agents Chemother.* 2007;51(8):2709–2715.
    <https://pubmed.ncbi.nlm.nih.gov/17517845/>
    — The model's tafenoquine PK (`VTQ` about 1600 L, half-life about 15 days).

40. **Commons RJ, Simpson JA, Thriemer K, et al.** The haematological
    consequences of *Plasmodium vivax* malaria after chloroquine treatment
    with and without primaquine: a WorldWide Antimalarial Resistance Network
    systematic review and individual patient data meta-analysis.
    *BMC Med.* 2019;17(1):151.
    <https://pubmed.ncbi.nlm.nih.gov/31366382/>

---

## 7. Methaemoglobin · methylene blue — the second consumer of the same NADPH pool

**In the model, the contraindication to methylene blue is not a fact to be memorised but the outcome of a calculation.**

41. **Rosen PJ, Johnson C, McGehee WG, Beutler E.** Failure of methylene blue
    treatment in toxic methemoglobinemia: association with glucose-6-phosphate
    dehydrogenase deficiency. *Ann Intern Med.* 1971;75(1):83–86.
    <https://pubmed.ncbi.nlm.nih.gov/5091568/>
    — The original clinical observation of scenario 7.

42. **Wright RO, Lewander WJ, Woolf AD.** Methemoglobinemia: etiology,
    pharmacology, and clinical management.
    *Ann Emerg Med.* 1999;34(5):646–656.
    <https://pubmed.ncbi.nlm.nih.gov/10533013/>

43. **Cortazzo JA, Lichtman AD.** Methemoglobinemia: a review and
    recommendations for management.
    *J Cardiothorac Vasc Anesth.* 2014;28(4):1043–1047.
    <https://pubmed.ncbi.nlm.nih.gov/24075607/>

44. **Coleman MD, Coleman NA.** Drug-induced methaemoglobinaemia: treatment
    issues. *Drug Saf.* 1996;14(6):394–405.
    <https://pubmed.ncbi.nlm.nih.gov/8828017/>

45. **Coleman MD, Jacobus DP.** Reduction of dapsone hydroxylamine to dapsone
    during methaemoglobin formation in human erythrocytes in vitro.
    *Biochem Pharmacol.* 1993;45(5):1027–1033.
    <https://pubmed.ncbi.nlm.nih.gov/8461032/>
    — The model's `KMETDAP`. The co-oxidation cycle in which one hydroxylamine molecule catalytically oxidises
    many haemoglobin molecules.

46. **Zuidema J, Hilbers-Modderman ES, Merkus FW.** Clinical pharmacokinetics
    of dapsone. *Clin Pharmacokinet.* 1986;11(4):299–315.
    <https://pubmed.ncbi.nlm.nih.gov/3530584/>
    — The model's dapsone PK (`VDP`, `CLDP`, half-life 20–30 hours).

---

## 8. Rasburicase · tumour lysis syndrome — the tumour sets the oxidant dose

47. **Relling MV, McDonagh EM, Chang T, et al.** Clinical Pharmacogenetics
    Implementation Consortium (CPIC) guidelines for rasburicase therapy in the
    context of G6PD deficiency genotype.
    *Clin Pharmacol Ther.* 2014;96(2):169–174.
    <https://pubmed.ncbi.nlm.nih.gov/24787449/>

48. **Browning LA, Kruse JA.** Hemolysis and methemoglobinemia secondary to
    rasburicase administration.
    *Ann Pharmacother.* 2005;39(11):1932–1935.
    <https://pubmed.ncbi.nlm.nih.gov/16204389/>
    — How the stoichiometry of one mole of hydrogen peroxide per mole of uric acid appears in the clinic.

49. **Sherwood GB, Paschal RD, Adamski J.** Rasburicase-induced
    methemoglobinemia: case report, literature review, and proposed treatment
    algorithm. *Clin Case Rep.* 2016;4(4):315–319.
    <https://pubmed.ncbi.nlm.nih.gov/27099721/>

50. **Coiffier B, Altman A, Pui C-H, Younes A, Cairo MS.** Guidelines for the
    management of pediatric and adult tumor lysis syndrome: an evidence-based
    review. *J Clin Oncol.* 2008;26(16):2767–2778.
    <https://pubmed.ncbi.nlm.nih.gov/18509186/>
    — The size of the uric acid load in scenario 8 (15–25 mg/dL).

51. **Cheuk DKL, Chiang AKS, Chan GCF, Ha SY.** Urate oxidase for the
    prevention and treatment of tumour lysis syndrome in children with cancer.
    *Cochrane Database Syst Rev.* 2017;3:CD006945.
    <https://pubmed.ncbi.nlm.nih.gov/28272834/>

---

## 9. Neonatal jaundice · kernicterus — it is multiplication, not addition

52. **Kaplan M, Renbaum P, Levy-Lahad E, Hammerman C, Lahad A, Beutler E.**
    Gilbert syndrome and glucose-6-phosphate dehydrogenase deficiency: a
    dose-dependent genetic interaction crucial to neonatal hyperbilirubinemia.
    *Proc Natl Acad Sci USA.* 1997;94(22):12128–12132.
    <https://pubmed.ncbi.nlm.nih.gov/9342374/>
    — **The source paper for scenario 10.** It showed directly the interaction in which either defect alone is safe
    and only the two together become dangerous. Why the model writes production and conjugation as
    numerator and denominator.

53. **Kaplan M, Hammerman C.** Glucose-6-phosphate dehydrogenase deficiency
    and severe neonatal hyperbilirubinemia: a complexity of interactions
    between genes and environment.
    *Semin Fetal Neonatal Med.* 2010;15(3):148–156.
    <https://pubmed.ncbi.nlm.nih.gov/19942984/>

54. **Watchko JF, Tiribelli C.** Bilirubin-induced neurologic damage —
    mechanisms and management approaches.
    *N Engl J Med.* 2013;369(21):2021–2030.
    <https://pubmed.ncbi.nlm.nih.gov/24256380/>
    — That free bilirubin (`BFREE`), not total bilirubin, is the species that crosses into the brain.
    Why the model solves albumin binding explicitly.

55. **Bhutani VK, Johnson LH.** Kernicterus in the 21st century: frequently
    asked questions. *J Perinatol.* 2009;29(Suppl 1):S20–S24.
    <https://pubmed.ncbi.nlm.nih.gov/19177056/>

56. **Slusher TM, Zamora TG, Appiah D, et al.** Burden of severe neonatal
    jaundice: a systematic review and meta-analysis.
    *BMJ Paediatr Open.* 2017;1(1):e000105.
    <https://pubmed.ncbi.nlm.nih.gov/29637134/>
    — The epidemiological evidence that G6PD deficiency is a leading cause of kernicterus worldwide.

57. **Bosma PJ, Chowdhury JR, Bakker C, et al.** The genetic basis of the
    reduced expression of bilirubin UDP-glucuronosyltransferase 1 in
    Gilbert's syndrome. *N Engl J Med.* 1995;333(18):1171–1175.
    <https://pubmed.ncbi.nlm.nih.gov/7565971/>
    — The (TA)n promoter genotype behind the model's `FUGT`.

---

## 10. Intravascular haemolysis · free haemoglobin · acute kidney injury

58. **Rother RP, Bell L, Hillmen P, Gladwin MT.** The clinical sequelae of
    intravascular hemolysis and extracellular plasma hemoglobin: a novel
    mechanism of human disease. *JAMA.* 2005;293(13):1653–1662.
    <https://pubmed.ncbi.nlm.nih.gov/15811985/>

59. **Schaer DJ, Buehler PW, Alayash AI, Belcher JD, Vercellotti GM.**
    Hemolysis and free hemoglobin revisited: exploring hemoglobin and hemin
    scavengers as a novel class of therapeutic proteins.
    *Blood.* 2013;121(8):1276–1284.
    <https://pubmed.ncbi.nlm.nih.gov/23264591/>
    — The model's haptoglobin saturation (`HP0`, `STOIHP`) and the renal filtration that follows it.

60. **Van Avondt K, Nur E, Zeerleder S.** Mechanisms of haemolysis-induced
    kidney injury. *Nat Rev Nephrol.* 2019;15(11):671–692.
    <https://pubmed.ncbi.nlm.nih.gov/31455889/>
    — The model's `TUB` → `GFRF` → `CREA` axis.

---

## 11. Red cell lifespan · erythropoiesis kinetics (the structural skeleton of the model)

61. **Franco RS.** Measurement of red cell lifespan and aging.
    *Transfus Med Hemother.* 2012;39(5):302–307.
    <https://pubmed.ncbi.nlm.nih.gov/23801920/>
    — The 120-day lifespan and its distribution. The basis for the model using a chain of 12 transit compartments (Erlang),
    and the limits of that approximation.

62. **Higgins JM, Mahadevan L.** Physiological and pathological population
    dynamics of circulating human red blood cells.
    *Proc Natl Acad Sci USA.* 2010;107(47):20587–20592.
    <https://pubmed.ncbi.nlm.nih.gov/21059904/>
    — A quantitative framework treating red cells as an age-structured population.

63. **Cohen RM, Franco RS, Khera PK, et al.** Red cell life span heterogeneity
    in hematologically normal people is sufficient to alter HbA1c.
    *Blood.* 2008;112(10):4284–4291.
    <https://pubmed.ncbi.nlm.nih.gov/18694998/>

64. **Krzyzanski W, Perez-Ruixo JJ.** An assessment of recombinant human
    erythropoietin effect on reticulocyte production rate and lifespan
    distribution in healthy subjects.
    *Pharm Res.* 2007;24(4):758–772.
    <https://pubmed.ncbi.nlm.nih.gov/17318416/>
    — The model's EPO → production rate · shortened marrow transit time · shift reticulocyte structure.

65. **Fisher JW.** Erythropoietin: physiology and pharmacology update.
    *Exp Biol Med (Maywood).* 2003;228(1):1–14.
    <https://pubmed.ncbi.nlm.nih.gov/12524467/>
    — The log-linear response of EPO to the haemoglobin deficit (`GEPO`).

---

## 12. Diagnosis — the test lies precisely when it is needed

66. **Domingo GJ, Satyagraha AW, Anvikar A, et al.** G6PD testing in support
    of treatment and elimination of malaria: recommendations for evaluation
    of G6PD tests. *Malar J.* 2013;12:391.
    <https://pubmed.ncbi.nlm.nih.gov/24188096/>

67. **LaRue N, Kahn M, Murray M, et al.** Comparison of quantitative and
    qualitative tests for glucose-6-phosphate dehydrogenase deficiency.
    *Am J Trop Med Hyg.* 2014;91(4):854–861.
    <https://pubmed.ncbi.nlm.nih.gov/25071003/>

68. **Bancone G, Chu CS, Chowwiwat N, et al.** Suitability of capillary blood
    for quantitative assessment of G6PD activity and performances of G6PD
    point-of-care tests. *Am J Trop Med Hyg.* 2015;92(4):818–824.
    <https://pubmed.ncbi.nlm.nih.gov/25646252/>

69. **Beutler E, Blume KG, Kaplan JC, Löhr GW, Ramot B, Valentine WN.**
    International Committee for Standardization in Haematology: recommended
    methods for red-cell enzyme analysis.
    *Br J Haematol.* 1977;35(2):331–340.
    <https://pubmed.ncbi.nlm.nih.gov/857853/>
    — The fact that the standard quantitative assay reports the **whole-blood mean activity**. The model's `G6PDPCT`
    mimics exactly this quantity, and shows how a rise in reticulocytes pulls it up into the normal
    range (scenario 9).

70. **Ley B, Winasti Satyagraha A, Rahmat H, et al.** Performance of the
    Access Bio/CareStart rapid diagnostic test for the detection of
    glucose-6-phosphate dehydrogenase deficiency: a systematic review and
    meta-analysis. *PLoS Med.* 2019;16(12):e1002992.
    <https://pubmed.ncbi.nlm.nih.gov/31834890/>

---

## 13. Drug safety lists · prescribing guidance

71. **Gammal RS, Pirmohamed M, Somogyi AA, et al.** Expanded Clinical
    Pharmacogenetics Implementation Consortium guideline for medication use
    in the context of G6PD genotype.
    *Clin Pharmacol Ther.* 2023;113(5):973–985.
    <https://pubmed.ncbi.nlm.nih.gov/36049896/>
    — The current CPIC prescribing guideline. The list the model's drug selection follows.

72. **Youngster I, Arcavi L, Schechmaster R, et al.** Medications and
    glucose-6-phosphate dehydrogenase deficiency: an evidence-based review.
    *Drug Saf.* 2010;33(9):713–726.
    <https://pubmed.ncbi.nlm.nih.gov/20701405/>
    — A rare paper addressing head-on the point that much of the "list of contraindicated drugs" rests on
    weak evidence.

73. **World Health Organization.** *WHO Guidelines for Malaria.* Geneva: WHO.
    <https://www.who.int/publications/i/item/guidelines-for-malaria>
    — The recommendation of primaquine 45 mg once weekly × 8 weeks, and the requirement for a quantitative assay before tafenoquine.

74. **US FDA.** KRINTAFEL (tafenoquine) tablets — prescribing information.
    <https://www.accessdata.fda.gov/drugsatfda_docs/label/2018/210795s000lbl.pdf>
    — The source document for the requirement of a **quantitative** G6PD activity assay before dosing (above 70% of normal).

---

## 14. The chemistry of favism

75. **Arese P, De Flora A.** Pathophysiology of hemolysis in glucose-6-
    phosphate dehydrogenase deficiency.
    *Semin Hematol.* 1990;27(1):1–40.
    <https://pubmed.ncbi.nlm.nih.gov/2405497/>

76. **Chevion M, Navok T, Glaser G, Mager J.** The chemistry of favism-
    inducing compounds: the properties of isouramil and divicine and their
    reaction with glutathione.
    *Eur J Biochem.* 1982;127(2):405–409.
    <https://pubmed.ncbi.nlm.nih.gov/7140768/>
    — The model's `SFV`. Divicine and isouramil react directly with GSH and redox-cycle.

77. **Mavelli I, Ciriolo MR, Rotilio G, De Sole P, Castorino M, Stabile A.**
    Favism: a hemolytic disease associated with increased superoxide
    dismutase and decreased glutathione peroxidase activities in red blood
    cells. *Eur J Biochem.* 1984;139(1):13–18.
    <https://pubmed.ncbi.nlm.nih.gov/6698003/>

---

## 15. Tools

- **mrgsolve** — R-based ODE/PK-PD simulation: <https://mrgsolve.org/>
- **Graphviz** — rendering the mechanistic map: <https://graphviz.org/>
- **Shiny** — the interactive dashboard: <https://shiny.posit.co/>
- **The G6PD Deficiency Association drug list**: <https://www.g6pd.org/>
- **PharmGKB — G6PD**: <https://www.pharmgkb.org/gene/PA28469>

---

## Honest notes on this model's limitations

1. **No fitting was done.** The parameters were fixed at the central values of the literature above, and
   only the oxidant scale factors (`SPQ` and the like) were adjusted so that scenario 1 (the Dern curve) and
   scenario 10 (neonatal bilirubin) fell inside the reported ranges. No formal estimation or validation
   against individual patient data was carried out.
2. **The age structure is a 12-compartment Erlang approximation.** It is more spread out than the real red cell
   lifespan distribution (CV about 29%), so a\* becomes a narrow band rather than a knife edge. Age differences finer
   than 10 days cannot be resolved with this model.
3. **Glutathione is collapsed to a quasi-steady state.** A minute-scale process has been solved algebraically
   for a day-scale disease, so nothing can be said about the dynamics over the minutes to tens of minutes right after exposure.
4. **A heterozygous female is run twice and mixed** (`mix_mosaic()`). This is exact for the red cell
   populations, but the plasma indices (bilirubin, free Hb) are recombined by mass weighting and so are an
   approximation.
5. **The Class I (CNSHA) phenotype is milder than reality.** The baseline anaemia the model produces is in the
   mild range, and it does not reproduce the chronic transfusion dependence of a real Class I patient.
6. For education, research and hypothesis generation; **it must not be used for clinical decision-making.**
