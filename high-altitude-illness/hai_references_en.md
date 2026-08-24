# High-Altitude Illness (AMS · HACE · HAPE) — References

The source literature behind the high-altitude illness QSP model. Each section
corresponds to the matching component of the model, and references used **directly** for
parameter calibration are marked **calibration**.

> **PMID verification.** Every PubMed ID below was looked up through the NCBI E-utilities (`esearch` +
> `esummary`) and its **title, authors, journal and year actually confirmed**.
> Not one PMID here was written from memory. Candidates whose lookup turned out to point at
> a paper other than the one intended were dropped from the list.

---

## 1. The physics of altitude — barometric pressure and inspired oxygen partial pressure

The model's only exogenous variable is `PIO2 = FiO2 × (PB(h) − 47)`. Everything else
follows from it.

1. West JB. **Prediction of barometric pressures at high altitude with the use
   of model atmospheres.** J Appl Physiol. 1996;81(4):1850.
   <https://pubmed.ncbi.nlm.nih.gov/8904608/> — **calibration** (the model's
   barometric-pressure equation `PB = exp(6.63268 − 0.1112h − 0.00149h²)`)
2. West JB. **Barometric pressures on Mt. Everest: new data and physiological
   significance.** J Appl Physiol. 1999;86(3):1062.
   <https://pubmed.ncbi.nlm.nih.gov/10066724/> — **calibration** (summit 253 mmHg)
3. West JB, Hackett PH, Maret KH, et al. **Maximal exercise at extreme altitudes
   on Mount Everest.** J Appl Physiol. 1983;55(3):688.
   <https://pubmed.ncbi.nlm.nih.gov/6415008/>
4. Grocott MPW, Martin DS, Levett DZH, et al. **Arterial blood gases and oxygen
   content in climbers on Mount Everest.** N Engl J Med. 2009;360(2):140.
   <https://pubmed.ncbi.nlm.nih.gov/19129527/> — **calibration** (8400 m
   PaO2 24.6 mmHg · PaCO2 13.3 mmHg · SaO2 54 %)
5. West JB. **High-altitude medicine.** Am J Respir Crit Care Med.
   2012;186(12):1229. <https://pubmed.ncbi.nlm.nih.gov/23103737/>
6. Grocott M, Montgomery H, Vercueil A. **High-altitude physiology and
   pathophysiology: implications and relevance for intensive care medicine.**
   Crit Care. 2007;11(1):203. <https://pubmed.ncbi.nlm.nih.gov/17291330/>

## 2. Pulmonary gas exchange — alveolar gas equation · diffusion · shunt

7. Severinghaus JW. **Simple, accurate equations for human blood O2 dissociation
   computations.** J Appl Physiol. 1979;46(3):599.
   <https://pubmed.ncbi.nlm.nih.gov/35496/> — **calibration** (oxygen dissociation curve
   `S = 1/(23400/(P³+150P)+1)`)
8. Piiper J, Scheid P. **Model for capillary-alveolar equilibration with special
   reference to O2 uptake in hypoxia.** Respir Physiol. 1981;46(3):193.
   <https://pubmed.ncbi.nlm.nih.gov/6798659/> — **calibration** (diffusion limitation
   `exp(−D_L/βQ̇)`; why the model's A–a gradient *narrows* from 9 mmHg at sea
   level to 3 mmHg at the summit)
9. Wagner PD, Gale GE, Moon RE, et al. **Pulmonary gas exchange in humans
   exercising at sea level and simulated altitude.** J Appl Physiol.
   1986;61(1):260. <https://pubmed.ncbi.nlm.nih.gov/3090012/>
10. Sutton JR, Reeves JT, Wagner PD, et al. **Operation Everest II: oxygen
    transport during exercise at extreme simulated altitude.** J Appl Physiol.
    1988;64(4):1309. <https://pubmed.ncbi.nlm.nih.gov/3132445/>
11. Calbet JAL, Boushel R, Rådegran G, et al. **Determinants of maximal oxygen
    uptake in severe acute hypoxia.** Am J Physiol Regul Integr Comp Physiol.
    2003;284(2):R291. <https://pubmed.ncbi.nlm.nih.gov/12388461/>
12. Calbet JAL, Lundby C. **Air to muscle O2 delivery during exercise at
    altitude.** High Alt Med Biol. 2009;10(2):123.
    <https://pubmed.ncbi.nlm.nih.gov/19555296/>
13. Wehrlin JP, Hallén J. **Linear decrease in V̇O2max and performance with
    increasing altitude in endurance athletes.** Eur J Appl Physiol.
    2006;96(4):404. <https://pubmed.ncbi.nlm.nih.gov/16311764/>

## 3. Ventilatory control and ventilatory acclimatisation

Defence 1 of the model. The key point is that the limit is set not by the chemoreceptors but by **acid-base**.

14. Dempsey JA, Forster HV. **Mediation of ventilatory adaptations.** Physiol
    Rev. 1982;62(1):262. <https://pubmed.ncbi.nlm.nih.gov/6798585/>
15. Dempsey JA, Smith CA. **Pathophysiology of human ventilatory control.**
    Eur Respir J. 2014;44(2):495.
    <https://pubmed.ncbi.nlm.nih.gov/24925922/>
16. Teppema LJ, Dahan A. **The ventilatory response to hypoxia in mammals:
    mechanisms, measurement, and analysis.** Physiol Rev. 2010;90(2):675.
    <https://pubmed.ncbi.nlm.nih.gov/20393196/> — **calibration** (HVR
    inter-individual range, multiplicative form of the peripheral–central interaction)
17. Weil JV, Byrne-Quinn E, Sodal IE, et al. **Hypoxic ventilatory drive in
    normal man.** J Clin Invest. 1970;49(6):1061.
    <https://pubmed.ncbi.nlm.nih.gov/5422012/> — **calibration**
    (isocapnic HVR slope → model `GP`)
18. Rahn H, Otis AB. **Man's respiratory response during and after
    acclimatization to high altitude.** Am J Physiol. 1949;157(3):445.
    <https://pubmed.ncbi.nlm.nih.gov/18151752/> — **calibration**
    (the classical PaO2–PaCO2 altitude line)
19. Guluzade NA, et al. **A test of the interaction between central and
    peripheral respiratory chemoreflexes in humans.** J Physiol.
    2023;601(19):4341. <https://pubmed.ncbi.nlm.nih.gov/37566804/>
20. Fan JL, Subudhi AW, Evero O, et al. **AltitudeOmics: enhanced
    cerebrovascular reactivity and ventilatory response to CO2 with
    high-altitude acclimatization and re-exposure.** J Appl Physiol.
    2014;116(7):911. <https://pubmed.ncbi.nlm.nih.gov/24356520/>
21. Forster HV, et al. **Control of breathing during exercise.** Compr Physiol.
    2012;2(1):743. <https://pubmed.ncbi.nlm.nih.gov/23728984/>

## 4. Acid-base — the rate-limiting step of acclimatisation

22. Krapf R, Beeler I, Hertner D, Hulter HN. **Chronic respiratory alkalosis.
    The effect of sustained hyperventilation on renal regulation of acid-base
    equilibrium.** N Engl J Med. 1991;324(20):1394.
    <https://pubmed.ncbi.nlm.nih.gov/1902283/> — **calibration** (chronic
    ΔHCO3⁻/ΔPaCO2 ≈ −0.5, renal time constant `TAU_REN`)
23. Siggaard-Andersen O. **Oxygen and acid-base parameters of arterial and mixed
    venous blood, relevant versus redundant.** Acta Anaesthesiol Scand Suppl.
    1995;107:21. <https://pubmed.ncbi.nlm.nih.gov/8599280/> — **calibration**
    (the Van Slyke equation → why the model reproduces acute buffering with no parameter of its own)

## 5. Acetazolamide — what it actually does

The model's conclusion: acetazolamide does not stimulate breathing. **It lowers the floor
(the apnoeic threshold) faster than it lowers the operating point.**

24. Swenson ER. **Carbonic anhydrase inhibitors and ventilation: a complex
    interplay of stimulation and suppression.** Eur Respir J. 1998;12(6):1242.
    <https://pubmed.ncbi.nlm.nih.gov/9877470/>
25. Swenson ER, Teppema LJ. **Prevention of acute mountain sickness by
    acetazolamide: as yet an unfinished story.** J Appl Physiol.
    2007;102(4):1305. <https://pubmed.ncbi.nlm.nih.gov/17194729/>
26. Leaf DE, Goldfarb DS. **Mechanisms of action of acetazolamide in the
    prophylaxis and treatment of acute mountain sickness.** J Appl Physiol.
    2007;102(4):1313. <https://pubmed.ncbi.nlm.nih.gov/17023566/>
27. Swenson ER, Hughes JMB. **Effects of acute and chronic acetazolamide on
    resting ventilation and ventilatory responses in men.** J Appl Physiol.
    1993;74(1):230. <https://pubmed.ncbi.nlm.nih.gov/8444696/> —
    **calibration** (magnitude of the HCO3⁻ fall → `ACZ_EMAX_REN`)
28. Teppema LJ, Balanos GM, Steinback CD, et al. **Effects of acetazolamide on
    ventilatory, cerebrovascular, and pulmonary vascular responses to hypoxia.**
    Am J Respir Crit Care Med. 2007;175(3):277.
    <https://pubmed.ncbi.nlm.nih.gov/17095745/>
29. Basnyat B, Gertsch JH, Johnson EW, et al. **Efficacy of low-dose
    acetazolamide (125 mg BID) for the prophylaxis of acute mountain sickness.**
    High Alt Med Biol. 2003;4(1):45.
    <https://pubmed.ncbi.nlm.nih.gov/12713711/> — **calibration** (the 125 mg dose arm)
30. Low EV, Avery AJ, Gupta V, et al. **Identifying the lowest effective dose of
    acetazolamide for the prophylaxis of acute mountain sickness: systematic
    review and meta-analysis.** BMJ. 2012;345:e6779.
    <https://pubmed.ncbi.nlm.nih.gov/23081689/> — **calibration** (RR)
31. Kayser B, Dumont L, Lysakowski C, et al. **Reappraisal of acetazolamide for
    the prevention of acute mountain sickness: a systematic review and
    meta-analysis.** High Alt Med Biol. 2012;13(2):82.
    <https://pubmed.ncbi.nlm.nih.gov/22724610/>
32. van Patot MCT, Leadbetter G, Keyes LE, et al. **Prophylactic low-dose
    acetazolamide reduces the incidence and severity of acute mountain
    sickness.** High Alt Med Biol. 2008;9(4):289.
    <https://pubmed.ncbi.nlm.nih.gov/19115912/>
33. Gertsch JH, Basnyat B, Johnson EW, Onopa J, Holck PS. **Randomised, double
    blind, placebo controlled comparison of ginkgo biloba and acetazolamide for
    prevention of acute mountain sickness.** BMJ. 2004;328(7443):797.
    <https://pubmed.ncbi.nlm.nih.gov/15070635/>
34. Richalet JP, Rivera M, Bouchet P, et al. **Acetazolamide: a treatment for
    chronic mountain sickness.** Am J Respir Crit Care Med. 2005;172(11):1427.
    <https://pubmed.ncbi.nlm.nih.gov/16126936/>

## 6. Periodic breathing during sleep — the CO2 reserve

The model's `CO2_reserve = PaCO2 − Bc` comes from here.

35. Dempsey JA. **Crossing the apnoeic threshold: causes and consequences.**
    Exp Physiol. 2005;90(1):13. <https://pubmed.ncbi.nlm.nih.gov/15572458/>
    — **calibration** (sea-level CO2 reserve 3–6 mmHg)
36. Dempsey JA, Veasey SC, Morgan BJ, O'Donnell CP. **Pathophysiology of sleep
    apnea.** Physiol Rev. 2010;90(1):47.
    <https://pubmed.ncbi.nlm.nih.gov/20086074/>
37. Khoo MCK, Kronauer RE, Strohl KP, Slutsky AS. **Factors inducing periodic
    breathing in humans: a general model.** J Appl Physiol. 1982;53(3):644.
    <https://pubmed.ncbi.nlm.nih.gov/7129986/> — (loop-gain formalism)
38. Bloch KE, Latshang TD, Turk AJ, et al. **Nocturnal periodic breathing during
    acclimatization at very high altitude at Mount Muztagh Ata (7,546 m).**
    Am J Respir Crit Care Med. 2010;182(4):562.
    <https://pubmed.ncbi.nlm.nih.gov/20442435/> — **calibration** (AHI versus altitude;
    ★ also the evidence for the nocturnal trend the model *fails to reproduce* — see the limitations section of the README)
39. Nussbaumer-Ochsner Y, Ursprung J, Siebenmann C, Maggiorini M, Bloch KE.
    **Effect of short-term acclimatization to high altitude on sleep and
    nocturnal breathing.** Sleep. 2012;35(3):419.
    <https://pubmed.ncbi.nlm.nih.gov/22379248/>
40. Fischer R, Lang SM, Leitl M, et al. **Theophylline and acetazolamide reduce
    sleep-disordered breathing at high altitude.** Eur Respir J. 2004;23(1):47.
    <https://pubmed.ncbi.nlm.nih.gov/14738230/> — **calibration**
    (magnitude of the AHI reduction with acetazolamide)
41. Burgess A, et al. **Loop gain response to increased cerebral blood flow at
    high altitude.** Sleep Breath. 2024;28(2):873.
    <https://pubmed.ncbi.nlm.nih.gov/38085496/>

## 7. Acute mountain sickness (AMS) — epidemiology · scoring · ascent rate

42. Roach RC, Hackett PH, Oelz O, et al. **The 2018 Lake Louise Acute Mountain
    Sickness Score.** High Alt Med Biol. 2018;19(1):4.
    <https://pubmed.ncbi.nlm.nih.gov/29583031/> — **calibration**
    (the four LLS subdomains, AMS definition: headache ≥1 and total score ≥3)
43. Hackett PH, Rennie D, Levine HD. **The incidence, importance, and
    prophylaxis of acute mountain sickness.** Lancet. 1976;2(7996):1149.
    <https://pubmed.ncbi.nlm.nih.gov/62991/>
44. Hackett PH, Roach RC. **High-altitude illness.** N Engl J Med.
    2001;345(2):107. <https://pubmed.ncbi.nlm.nih.gov/11450659/>
45. Bärtsch P, Swenson ER. **Acute high-altitude illnesses.** N Engl J Med.
    2013;369(17):1666. <https://pubmed.ncbi.nlm.nih.gov/24152275/>
46. Luks AM, Swenson ER, Bärtsch P. **Acute high-altitude sickness.** Eur Respir
    Rev. 2017;26(143):160096. <https://pubmed.ncbi.nlm.nih.gov/28143879/>
47. Maggiorini M, Bühler B, Walter M, Oelz O. **Prevalence of acute mountain
    sickness in the Swiss Alps.** BMJ. 1990;301(6756):853.
    <https://pubmed.ncbi.nlm.nih.gov/2282425/> — **calibration** (prevalence by altitude)
48. Schneider M, Bernasch D, Weymann J, Holle R, Bärtsch P. **Acute mountain
    sickness: influence of susceptibility, preexposure, and ascent rate.**
    Med Sci Sports Exerc. 2002;34(12):1886.
    <https://pubmed.ncbi.nlm.nih.gov/12471292/> — **calibration**
    (the comparator for the ascent-rate sweep)
49. Bloch KE, Turk AJ, Maggiorini M, et al. **Effect of ascent protocol on acute
    mountain sickness and success at Muztagh Ata, 7546 m.** High Alt Med Biol.
    2009;10(1):25. <https://pubmed.ncbi.nlm.nih.gov/19326598/>
50. Bärtsch P, Bailey DM, Berger MM, Knauth M, Baumgartner RW. **Acute mountain
    sickness: controversies and advances.** High Alt Med Biol. 2004;5(2):110.
    <https://pubmed.ncbi.nlm.nih.gov/15265333/>
51. Imray C, Wright A, Subudhi A, Roach R. **Acute mountain sickness:
    pathophysiology, prevention, and treatment.** Prog Cardiovasc Dis.
    2010;52(6):467. <https://pubmed.ncbi.nlm.nih.gov/20417340/>
52. Canoui-Poitrine F, Veerabudun K, Larmignat P, et al. **Risk prediction score
    for severe high altitude illness: a cohort study.** PLoS One.
    2014;9(7):e100642. <https://pubmed.ncbi.nlm.nih.gov/25068815/>
53. Luks AM, Auerbach PS, Freer L, et al. **Wilderness Medical Society Clinical
    Practice Guidelines for the Prevention and Treatment of Acute Altitude
    Illness: 2019 Update.** Wilderness Environ Med. 2019;30(4S):S3.
    <https://pubmed.ncbi.nlm.nih.gov/31248818/> — **calibration**
    (recommended ascent rates, prophylactic drug doses)

## 8. The brain — cerebral blood flow · oedema · intracranial pressure

Defence 3 of the model, and the quantitative form of the "tight-fit" hypothesis.

54. Wilson MH, Newman S, Imray CH. **The cerebral effects of ascent to high
    altitudes.** Lancet Neurol. 2009;8(2):175.
    <https://pubmed.ncbi.nlm.nih.gov/19161909/>
55. Ross RT. **The random nature of cerebral mountain sickness.** Lancet.
    1985;1(8435):990. <https://pubmed.ncbi.nlm.nih.gov/2859454/> —
    (the original "tight-fit" hypothesis → inter-individual variation in the model's `PVI`)
56. Marmarou A, Shulman K, Rosende RM. **A nonlinear analysis of the
    cerebrospinal fluid system and intracranial pressure dynamics.**
    J Neurosurg. 1978;48(3):332. <https://pubmed.ncbi.nlm.nih.gov/632857/> —
    **calibration** (pressure-volume index PVI, `ICP = ICP₀·10^(ΔV/PVI)`)
57. Lawley JS, Levine BD, Williams MA, et al. **Cerebral spinal fluid dynamics:
    effect of hypoxia and implications for high-altitude illness.** J Appl
    Physiol. 2016;120(2):251. <https://pubmed.ncbi.nlm.nih.gov/26494441/>
58. Hackett PH, Yarnell PR, Hill R, et al. **High-altitude cerebral edema
    evaluated with magnetic resonance imaging: clinical correlation and
    pathophysiology.** JAMA. 1998;280(22):1920.
    <https://pubmed.ncbi.nlm.nih.gov/9851477/> — **calibration**
    (vasogenic oedema, splenium of the corpus callosum)
59. Kallenberg K, Bailey DM, Christ S, et al. **Magnetic resonance imaging
    evidence of cytotoxic cerebral edema in acute mountain sickness.**
    J Cereb Blood Flow Metab. 2007;27(5):1064.
    <https://pubmed.ncbi.nlm.nih.gov/17024110/>
60. Schoonman GG, Sándor PS, Nirkko AC, et al. **Hypoxia-induced acute mountain
    sickness is associated with intracellular cerebral edema: a 3 T MRI study.**
    J Cereb Blood Flow Metab. 2008;28(1):198.
    <https://pubmed.ncbi.nlm.nih.gov/17519973/>
61. Sagoo RS, Hutchinson CE, Wright A, et al. **Magnetic resonance investigation
    into the mechanisms involved in the development of high-altitude cerebral
    edema.** J Cereb Blood Flow Metab. 2017;37(1):319.
    <https://pubmed.ncbi.nlm.nih.gov/26746867/>
62. Bailey DM, Bärtsch P, Knauth M, Baumgartner RW. **Emerging concepts in acute
    mountain sickness and high-altitude cerebral edema: from the molecular to
    the morphological.** Cell Mol Life Sci. 2009;66(22):3583.
    <https://pubmed.ncbi.nlm.nih.gov/19763397/>
63. Van Osta A, Moraine JJ, Mélot C, et al. **Effects of high altitude exposure
    on cerebral hemodynamics in normal subjects.** Stroke. 2005;36(3):557.
    <https://pubmed.ncbi.nlm.nih.gov/15692117/> — **calibration** (magnitude of the acute CBF rise)
64. Ainslie PN, Subudhi AW. **Cerebral blood flow at high altitude.** High Alt
    Med Biol. 2014;15(2):133. <https://pubmed.ncbi.nlm.nih.gov/24971767/>
65. Willie CK, Macleod DB, Shaw AD, et al. **Regional brain blood flow in man
    during acute changes in arterial blood gases.** J Physiol. 2012;590(14):3261.
    <https://pubmed.ncbi.nlm.nih.gov/22495584/> — **calibration**
    (the **saturating** form of CO2 reactivity → why the model does not use a linear reactivity)
66. Ainslie PN, Ogoh S, Burgess K, et al. **Differential effects of acute
    hypoxia and high altitude on cerebral blood flow velocity and dynamic
    cerebral autoregulation.** J Appl Physiol. 2008;104(2):490.
    <https://pubmed.ncbi.nlm.nih.gov/18048592/>
67. Subudhi AW, Fan JL, Evero O, et al. **AltitudeOmics: effect of ascent and
    acclimatization to 5260 m on regional cerebral oxygen delivery.** Exp
    Physiol. 2014;99(5):772. <https://pubmed.ncbi.nlm.nih.gov/24243839/>
68. Wilson MH, Imray CHE, Hargens AR. **The headache of high altitude and
    microgravity — similarities with clinical syndromes of cerebral venous
    hypertension.** High Alt Med Biol. 2011;12(4):379.
    <https://pubmed.ncbi.nlm.nih.gov/22087727/>
69. Wilson MH, Davagnanam I, Holland G, et al. **Cerebral venous system and
    anatomical predisposition to high-altitude headache.** Ann Neurol.
    2013;73(3):381. <https://pubmed.ncbi.nlm.nih.gov/23444324/>

## 9. High-altitude pulmonary oedema (HAPE) — capillary pressure and stress failure

How Defence 2 of the model becomes a disease. The key reference is number 62 (Maggiorini 2001):
HAPE begins not in inflammation but in **pressure**.

70. Maggiorini M, Mélot C, Pierre S, et al. **High-altitude pulmonary edema is
    initially caused by an increase in capillary pressure.** Circulation.
    2001;103(16):2078. <https://pubmed.ncbi.nlm.nih.gov/11319198/> —
    **calibration** (HAPE capillary pressure 19–26 mmHg versus ~13 mmHg in controls →
    `PCAP_CRIT`)
71. West JB, Mathieu-Costello O. **Stress failure of pulmonary capillaries: role
    in lung and heart disease.** Lancet. 1992;340(8822):762.
    <https://pubmed.ncbi.nlm.nih.gov/1356184/> — **calibration**
    (the supralinear pressure dependence of stress failure → `J ∝ (Pcap − Pcrit)^1.5`)
72. Elliott AR, Fu Z, Tsukimoto K, Prediletto R, Mathieu-Costello O, West JB.
    **Short-term reversibility of ultrastructural changes in pulmonary
    capillaries caused by stress failure.** J Appl Physiol. 1992;73(3):1150.
    <https://pubmed.ncbi.nlm.nih.gov/1400030/>
73. Bärtsch P, Mairbäurl H, Maggiorini M, Swenson ER. **Physiological aspects of
    high-altitude pulmonary edema.** J Appl Physiol. 2005;98(3):1101.
    <https://pubmed.ncbi.nlm.nih.gov/15703168/>
74. Swenson ER, Bärtsch P. **High-altitude pulmonary edema.** Compr Physiol.
    2012;2(4):2753. <https://pubmed.ncbi.nlm.nih.gov/23720264/>
75. Swenson ER, et al. **Early hours in the development of high-altitude
    pulmonary edema: time course and mechanisms.** J Appl Physiol.
    2020;128(6):1539. <https://pubmed.ncbi.nlm.nih.gov/32213112/> —
    (evidence that inflammation is a consequence and not a cause)
76. Dehnert C, Berger MM, Mairbäurl H, Bärtsch P. **High altitude pulmonary
    edema: a pressure-induced leak.** Respir Physiol Neurobiol.
    2007;158(2-3):266. <https://pubmed.ncbi.nlm.nih.gov/17602898/>
77. Hopkins SR, Garg J, Bolar DS, Balouch J, Levin DL. **Pulmonary blood flow
    heterogeneity during hypoxia and high-altitude pulmonary edema.** Am J
    Respir Crit Care Med. 2005;171(1):83.
    <https://pubmed.ncbi.nlm.nih.gov/15486339/> — **calibration**
    (★ direct evidence for the model's key parameter `a` = HPV heterogeneity)
78. Hultgren HN. **High-altitude pulmonary edema: current concepts.** Annu Rev
    Med. 1996;47:267. <https://pubmed.ncbi.nlm.nih.gov/8712781/>
79. Maggiorini M. **Prevention and treatment of high-altitude pulmonary edema.**
    Prog Cardiovasc Dis. 2010;52(6):500.
    <https://pubmed.ncbi.nlm.nih.gov/20417343/>
80. Grünig E, Mereles D, Hildebrandt W, et al. **Stress Doppler
    echocardiography for identification of susceptibility to high altitude
    pulmonary edema.** J Am Coll Cardiol. 2000;35(4):980.
    <https://pubmed.ncbi.nlm.nih.gov/10732898/> — **calibration**
    (exaggerated HPV in HAPE-susceptible subjects → `LAMBDA_MAX_S`)
81. Sartori C, Allemann Y, Trueb L, et al. **Augmented vasoreactivity in adult
    life associated with perinatal vascular insult.** Lancet.
    1999;353(9171):2205. <https://pubmed.ncbi.nlm.nih.gov/10392986/>
82. Allemann Y, Hutter D, Lipp E, et al. **Patent foramen ovale and
    high-altitude pulmonary edema.** JAMA. 2006;296(24):2954.
    <https://pubmed.ncbi.nlm.nih.gov/17190896/>
83. Duplain H, Sartori C, Lepori M, et al. **Exhaled nitric oxide in
    high-altitude pulmonary edema: role in the regulation of pulmonary vascular
    tone.** Am J Respir Crit Care Med. 2000;162(1):221.
    <https://pubmed.ncbi.nlm.nih.gov/10903245/>
84. Berger MM, Dehnert C, Bailey DM, et al. **Transpulmonary plasma ET-1 and
    nitrite differences in high altitude pulmonary hypertension.** High Alt Med
    Biol. 2009;10(1):17. <https://pubmed.ncbi.nlm.nih.gov/19326597/>

## 10. Hypoxic pulmonary vasoconstriction

85. Sylvester JT, Shimoda LA, Aaronson PI, Ward JPT. **Hypoxic pulmonary
    vasoconstriction.** Physiol Rev. 2012;92(1):367.
    <https://pubmed.ncbi.nlm.nih.gov/22298659/> — **calibration**
    (the PAO2-dependence curve of HPV → `P50_HPV`, `N_HPV`)
86. Weir EK, López-Barneo J, Buckler KJ, Archer SL. **Acute oxygen-sensing
    mechanisms.** N Engl J Med. 2005;353(19):2042.
    <https://pubmed.ncbi.nlm.nih.gov/16282179/>
87. Dorrington KL, Clar C, Young JD, et al. **Time course of the human pulmonary
    vascular response to 8 hours of isocapnic hypoxia.** Am J Physiol.
    1997;273(3):H1126. <https://pubmed.ncbi.nlm.nih.gov/9321798/> —
    **calibration** (the slow HPV component → `TAU_HPV_S`, `HPVS_MAX`)
88. Naeije R, Dedobbeleer C. **Pulmonary hypertension and the right ventricle in
    hypoxia.** Exp Physiol. 2013;98(8):1247.
    <https://pubmed.ncbi.nlm.nih.gov/23625956/>
89. Faoro V, et al. **Pulmonary circulation and gas exchange at exercise in
    Sherpas at high altitude.** J Appl Physiol. 2014;116(7):919.
    <https://pubmed.ncbi.nlm.nih.gov/23869067/>

## 11. Alveolar fluid clearance — the second positive feedback, in which the drain closes

90. Mairbäurl H. **Role of alveolar epithelial sodium transport in high altitude
    pulmonary edema (HAPE).** Respir Physiol Neurobiol. 2006;151(2-3):178.
    <https://pubmed.ncbi.nlm.nih.gov/16337225/>
91. Vivona ML, Matthay M, Chabaud MB, Friedlander G, Clerici C. **Hypoxia
    reduces alveolar epithelial sodium and fluid transport in rats: reversal by
    beta-adrenergic agonist treatment.** Am J Respir Cell Mol Biol.
    2001;25(5):554. <https://pubmed.ncbi.nlm.nih.gov/11713096/> —
    **calibration** (★ the model's second positive feedback `AFC_HYP50`)
92. Sartori C, Allemann Y, Duplain H, et al. **Salmeterol for the prevention of
    high-altitude pulmonary edema.** N Engl J Med. 2002;346(21):1631.
    <https://pubmed.ncbi.nlm.nih.gov/12023995/> — **calibration**
    (rise in clearance with a β2 agonist → `SAL_EMAX_AFC`)

## 12. HAPE prophylaxis and treatment trials

93. Bärtsch P, Maggiorini M, Ritter M, Noti C, Vock P, Oelz O. **Prevention of
    high-altitude pulmonary edema by nifedipine.** N Engl J Med.
    1991;325(18):1284. <https://pubmed.ncbi.nlm.nih.gov/1922223/> —
    **calibration** (HAPE incidence 63 % → 10 %)
94. Oelz O, Maggiorini M, Ritter M, et al. **Prevention and treatment of high
    altitude pulmonary edema by a calcium channel blocker.** Int J Sports Med.
    1992;13 Suppl 1:S65. <https://pubmed.ncbi.nlm.nih.gov/1483797/>
95. Maggiorini M, Brunner-La Rocca HP, Peth S, et al. **Both tadalafil and
    dexamethasone may reduce the incidence of high-altitude pulmonary edema: a
    randomized trial.** Ann Intern Med. 2006;145(7):497.
    <https://pubmed.ncbi.nlm.nih.gov/17015867/> — **calibration**
    (placebo 74 % · tadalafil 14 % · dexamethasone 29 % → the comparator for model
    scenarios 09–11)
96. Scherrer U, Vollenweider L, Delabays A, et al. **Inhaled nitric oxide for
    high-altitude pulmonary edema.** N Engl J Med. 1996;334(10):624.
    <https://pubmed.ncbi.nlm.nih.gov/8592525/>
97. Nieto Estrada VH, Molano Franco D, Medina RD, et al. **Interventions for
    preventing high altitude illness: Part 1. Commonly-used classes of drugs.**
    Cochrane Database Syst Rev. 2017;6:CD009761.
    <https://pubmed.ncbi.nlm.nih.gov/28653390/>

## 13. Dexamethasone · NSAIDs · physical measures

98. Ferrazzini G, Maggiorini M, Kriemler S, Bärtsch P, Oelz O. **Successful
    treatment of acute mountain sickness with dexamethasone.** Br Med J.
    1987;294(6584):1380. <https://pubmed.ncbi.nlm.nih.gov/3109663/>
99. Levine BD, Yoshimura K, Kobayashi T, et al. **Dexamethasone in the treatment
    of acute mountain sickness.** N Engl J Med. 1989;321(25):1707.
    <https://pubmed.ncbi.nlm.nih.gov/2687688/> — **calibration**
    (★ the primary data for the asymmetry in which symptoms improve while SaO2 does not move)
100. Burns P, et al. **Altitude sickness prevention with ibuprofen relative to
     acetazolamide.** Am J Med. 2019;132(2):247.
     <https://pubmed.ncbi.nlm.nih.gov/30419226/>
101. Wang J, et al. **Comparative effects of pharmacological interventions for
     the prevention of acute mountain sickness: a systematic review and network
     meta-analysis.** Travel Med Infect Dis. 2025;66:102850.
     <https://pubmed.ncbi.nlm.nih.gov/40383249/>
102. Bärtsch P, Merki B, Hofstetter D, Maggiorini M, Kayser B, Oelz O.
     **Treatment of acute mountain sickness by simulated descent: a randomised
     controlled trial.** BMJ. 1993;306(6885):1098.
     <https://pubmed.ncbi.nlm.nih.gov/8495155/> — **calibration**
     (Gamow bag = simulated descent; the clinical counterpart of the model's "descent-equivalent metres" calculation)
103. Freeman K, Shalit M, Stroh G. **Use of the Gamow Bag by EMT-basic park
     rangers for treatment of high-altitude pulmonary edema and high-altitude
     cerebral edema.** Wilderness Environ Med. 2004;15(3):198.
     <https://pubmed.ncbi.nlm.nih.gov/15473460/>
104. Luks AM, Swenson ER. **Travel to high altitude with pre-existing lung
     disease.** Eur Respir J. 2007;29(4):770.
     <https://pubmed.ncbi.nlm.nih.gov/17400877/>

## 14. Erythropoiesis · plasma volume · viscosity

The model's conclusion: the rise in Hct in week 1 is not red cells but **water**. And the optimal
haematocrit is 43 %, independent of altitude.

105. Eckardt KU, Boutellier U, Kurtz A, et al. **Rate of erythropoietin
     formation in humans in response to acute hypobaric hypoxia.** J Appl
     Physiol. 1989;66(4):1785. <https://pubmed.ncbi.nlm.nih.gov/2732171/> —
     **calibration** (EPO peak at 24–48 h, ×3–5)
106. Ge RL, Witkowski S, Zhang Y, et al. **Determinants of erythropoietin
     release in response to short-term hypobaric hypoxia.** J Appl Physiol.
     2002;92(6):2361. <https://pubmed.ncbi.nlm.nih.gov/12015348/>
107. Siebenmann C, Robach P, Lundby C. **Regulation of blood volume in
     lowlanders exposed to high altitude.** J Appl Physiol. 2017;123(4):957.
     <https://pubmed.ncbi.nlm.nih.gov/28572493/> — **calibration**
     (extent of plasma volume contraction → `PV_CONTRACT`)
108. Rasmussen P, Siebenmann C, Díaz V, Lundby C. **Red cell volume expansion at
     altitude: a meta-analysis and Monte Carlo simulation.** Med Sci Sports
     Exerc. 2013;45(9):1767. <https://pubmed.ncbi.nlm.nih.gov/23502972/> —
     **calibration** (rate of Hb mass increase ≈ +1 %/week)
109. Winslow RM, Monge CC, Statham NJ, et al. **Variability of oxygen affinity
     of blood: human subjects native to high altitude.** J Appl Physiol.
     1981;51(6):1411. <https://pubmed.ncbi.nlm.nih.gov/7319874/>
110. Mairbäurl H. **Red blood cells in sports: effects of exercise and training
     on oxygen supply by red blood cells.** Front Physiol. 2013;4:332.
     <https://pubmed.ncbi.nlm.nih.gov/24273518/>

## 15. Chronic mountain sickness and population adaptation

111. Villafuerte FC, Corante N. **Chronic mountain sickness: clinical aspects,
     etiology, management, and treatment.** High Alt Med Biol. 2016;17(2):61.
     <https://pubmed.ncbi.nlm.nih.gov/27218284/> — **calibration**
     (Qinghai score, diagnostic Hb threshold)
112. León-Velarde F, Maggiorini M, Reeves JT, et al. **Consensus statement on
     chronic and subacute high altitude diseases.** High Alt Med Biol.
     2005;6(2):147. <https://pubmed.ncbi.nlm.nih.gov/16060849/>
113. Beall CM. **Two routes to functional adaptation: Tibetan and Andean
     high-altitude natives.** Proc Natl Acad Sci USA. 2007;104 Suppl 1:8655.
     <https://pubmed.ncbi.nlm.nih.gov/17494744/> — (★ the anthropological counterpart
     of what the model's optimal haematocrit result is saying)
114. Beall CM. **Natural selection on EPAS1 (HIF2α) associated with low
     hemoglobin concentration in Tibetan highlanders.** Proc Natl Acad Sci USA.
     2010;107(25):11459. <https://pubmed.ncbi.nlm.nih.gov/20534544/>
115. Huerta-Sánchez E, Jin X, Asan, et al. **Altitude adaptation in Tibetans
     caused by introgression of Denisovan-like DNA.** Nature. 2014;512(7513):194.
     <https://pubmed.ncbi.nlm.nih.gov/25043035/>
116. Simonson TS, Yang Y, Huff CD, et al. **Genetic evidence for high-altitude
     adaptation in Tibet.** Science. 2010;329(5987):72.
     <https://pubmed.ncbi.nlm.nih.gov/20466884/>
117. Beall CM. **Human adaptability studies at high altitude: research designs
     and major concepts during fifty years of discovery.** Am J Hum Biol.
     2013;25(2):141. <https://pubmed.ncbi.nlm.nih.gov/23349118/>
118. Levett DZH, Radford EJ, Menassa DA, et al. **Acclimatization of skeletal
     muscle mitochondria to high-altitude hypoxia during an ascent of Everest.**
     FASEB J. 2012;26(4):1431. <https://pubmed.ncbi.nlm.nih.gov/22186874/>

## 16. Fluid and neurohumoral responses — the antidiuretic phenotype of AMS

119. Bärtsch P, Maggiorini M, Schobersberger W, et al. **Enhanced
     exercise-induced rise of aldosterone and vasopressin preceding mountain
     sickness.** J Appl Physiol. 1991;71(1):136.
     <https://pubmed.ncbi.nlm.nih.gov/1917735/> — **calibration**
     (AMS = antidiuretic phenotype → the model's `ADH`–`FLUID` pathway)
120. Loeppky JA, Icenogle MV, Maes D, et al. **Early fluid retention and severe
     acute mountain sickness.** J Appl Physiol. 2005;98(2):591.
     <https://pubmed.ncbi.nlm.nih.gov/15501929/>
121. Swenson ER, Duncan TB, Goldberg SV, Ramirez G, Ahmad S, Schoene RB.
     **Diuretic effect of acute hypoxia in humans: relationship to hypoxic
     ventilatory responsiveness and renal hormones.** J Appl Physiol.
     1995;78(2):377. <https://pubmed.ncbi.nlm.nih.gov/7759405/> —
     (the primary data for "those who acclimatise well pass urine")

---

## Component-to-reference map

| Model component | Equation/parameter | Reference |
|---|---|---|
| Barometric pressure | `PB = exp(6.63268 − 0.1112h − 0.00149h²)` | 1, 2 |
| Alveolar gas equation | `PAO2 = PIO2 − PaCO2·[FiO2+(1−FiO2)/R]` | 5, 18 |
| Oxygen dissociation curve | Severinghaus + Bohr + 2,3-DPG | 7, 109 |
| Diffusion limitation | Piiper–Scheid `exp(−D_L/βQ̇)` | 8, 9, 10 |
| Central apnoeic threshold `Bc` | `[HCO3⁻]_CSF/(0.03·10^(pH*−6.1))` | 35, 36 |
| Peripheral chemoreflex gain `GP` | multiplicative O₂×CO₂ interaction | 16, 17, 19 |
| Ventilatory acclimatisation `VAH` | carotid body plasticity, τ ≈ 60 h | 14, 16, 20 |
| Renal acid-base `TAU_REN` | τ ≈ 34 h, ΔHCO3⁻/ΔPaCO2 ≈ −0.5 | 22, 23 |
| Acetazolamide PD | renal CA + choroid plexus CA | 24–28 |
| CO₂ reserve | `PaCO2 − Bc` | 35, 38, 40 |
| HPV `λ(PAO2)` | Hill, P50 45 mmHg | 85, 86, 87 |
| HPV heterogeneity `a` | two-compartment split, upper bound = 1/(1−a) | 77, 80 |
| Capillary stress failure | `J ∝ (Pcap − 19.5)^1.5` | 70, 71, 72 |
| Alveolar fluid clearance `AFC` | hypoxia inhibits ENaC/Na-K-ATPase | 90, 91, 92 |
| Cerebral CO₂ reactivity | **saturating** sigmoid | 65, 66 |
| Monro–Kellie / PVI | `ICP = ICP₀·10^(ΔV/PVI)` | 55, 56, 57 |
| Vasogenic/cytotoxic oedema | VEGF·BBB / Na-K-ATPase failure | 58, 59, 60, 61 |
| Lake Louise score | headache · gastrointestinal · fatigue · dizziness (0–3 each) | 42 |
| Plasma volume contraction | −20 %, τ ≈ 40 h | 107 |
| Erythropoiesis | EPO peak at 24–48 h, Hb mass +1 %/week | 105, 106, 108 |
| Optimal Hct | `Hct* = 1/(k·γ)`, k = 2.31 | 108, 110, 113 |
| Nifedipine / tadalafil / dexamethasone | HAPE prophylaxis, 3-arm | 93, 95 |
| Salmeterol | clearance raised on its own | 92 |
| The dexamethasone asymmetry | symptoms↓, PaO2 unchanged | 99 |
| Gamow bag | +105 mmHg = simulated descent | 102, 103 |
| Chronic mountain sickness | Qinghai score, Hb threshold | 111, 112 |

---

## Where this model diverges from the literature

Set down honestly. The detailed numbers are in the "Limitations" section of `README_en.md`.

1. **Time course of nocturnal periodic breathing (reference 38).** The model explains
   instability with the CO₂ reserve alone, so it predicts that **the first night after arrival is the worst**
   and that things improve thereafter. Bloch et al.'s field polysomnography reported that at 4559 m
   the AHI was *higher* on day 3 than on day 1. The missing mechanism is the peripheral chemoreflex
   gain (loop gain) that rises with acclimatisation, which the model carries into ventilation
   but not into the stability index.
2. **PaCO2 at 8400 m (reference 4).** Model 14.6 mmHg versus observed 13.3 mmHg — at extreme
   altitude it underpredicts ventilation by about 10 %.
3. **Optimal haematocrit (reference 113).** The model yields 43 % under the assumption (γ = 1)
   that the circulation absorbs none of the viscosity burden. Real cardiac output is
   regulated, so this value should be read as a lower bound, and a γ sensitivity analysis is
   presented alongside it in the README.
