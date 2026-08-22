# Spontaneous Intracerebral Haemorrhage (ICH) — references

A list of the literature used for the structure, parameters and validation anchors of the QSP model
(`ich_qsp_model.dot`, `ich_mrgsolve_model.R`, `ich_shiny_app.R`).

**Every PMID was actually looked up through the PubMed E-utilities (esearch + esummary) and its title, authors,
journal and year cross-checked.** No PMID was written from memory, and any reference not confirmed by look-up
was excluded from the list.

The **▶ Model link** item at the end of each section states which equation, parameter or structural choice of the model
that group of references supports.

---

## 1. Overview · Epidemiology · Guidelines

1. Sheth KN. **Spontaneous Intracerebral Hemorrhage.** *N Engl J Med.* 2022;387:1589-1596.
   <https://pubmed.ncbi.nlm.nih.gov/36300975/>
2. Puy L, Parry-Jones AR, Sandset EC, Dowlatshahi D, Ziai W, Cordonnier C.
   **Intracerebral haemorrhage.** *Nat Rev Dis Primers.* 2023;9:14.
   <https://pubmed.ncbi.nlm.nih.gov/36928219/>
3. Greenberg SM, Ziai WC, Cordonnier C, et al. **2022 Guideline for the Management of
   Patients With Spontaneous Intracerebral Hemorrhage: A Guideline From the American
   Heart Association/American Stroke Association.** *Stroke.* 2022;53:e282-e361.
   <https://pubmed.ncbi.nlm.nih.gov/35579034/>
4. Cordonnier C, Demchuk A, Ziai W, Anderson CS. **Intracerebral haemorrhage: current
   approaches to acute management.** *Lancet.* 2018;392:1257-1268.
   <https://pubmed.ncbi.nlm.nih.gov/30319113/>
5. Magid-Bernstein J, Girard R, Polster S, et al. **Cerebral Hemorrhage: Pathophysiology,
   Treatment, and Future Directions.** *Circ Res.* 2022;130:1204-1229.
   <https://pubmed.ncbi.nlm.nih.gov/35420918/>
6. Keep RF, Hua Y, Xi G. **Intracerebral haemorrhage: mechanisms of injury and therapeutic
   targets.** *Lancet Neurol.* 2012;11:720-731.
   <https://pubmed.ncbi.nlm.nih.gov/22698888/>
7. Xi G, Keep RF, Hoff JT. **Mechanisms of brain injury after intracerebral haemorrhage.**
   *Lancet Neurol.* 2006;5:53-63. <https://pubmed.ncbi.nlm.nih.gov/16361023/>
8. Aronowski J, Zhao X. **Molecular pathophysiology of cerebral hemorrhage: secondary
   brain injury.** *Stroke.* 2011;42:1781-1786.
   <https://pubmed.ncbi.nlm.nih.gov/21527759/>
9. van Asch CJ, Luitse MJ, Rinkel GJ, van der Tweel I, Algra A, Klijn CJ. **Incidence,
   case fatality, and functional outcome of intracerebral haemorrhage over time,
   according to age, sex, and ethnic origin: a systematic review and meta-analysis.**
   *Lancet Neurol.* 2010;9:167-176. <https://pubmed.ncbi.nlm.nih.gov/20056489/>
10. GBD 2019 Stroke Collaborators. **Global, regional, and national burden of stroke and
    its risk factors, 1990-2019.** *Lancet Neurol.* 2021;20:795-820.
    <https://pubmed.ncbi.nlm.nih.gov/34487721/>
11. Qureshi AI, Tuhrim S, Broderick JP, Batjer HH, Hondo H, Hanley DF. **Spontaneous
    intracerebral hemorrhage.** *N Engl J Med.* 2001;344:1450-1460.
    <https://pubmed.ncbi.nlm.nih.gov/11346811/>
12. Meretoja A, Strbian D, Putaala J, et al. **SMASH-U: a proposal for etiologic
    classification of intracerebral hemorrhage.** *Stroke.* 2012;43:2592-2597.
    <https://pubmed.ncbi.nlm.nih.gov/22858729/>
13. Hostettler IC, Seiffge DJ, Werring DJ. **Intracerebral hemorrhage: an update on
    diagnosis and treatment.** *Expert Rev Neurother.* 2019;19:679-694.
    <https://pubmed.ncbi.nlm.nih.gov/31188036/>

**▶ Model link** — the skeleton of the overall structure, the deep/lobar distinction of cluster 1 (aetiology · small-vessel lesions),
the `LOC` parameter, the baseline incidence of the mortality logit (`MORTA`), and the design that takes 90-day mRS as
the primary outcome. Refs 6-8 are the basis of the "primary injury (mass effect) + secondary injury (haem · iron · inflammation)" dichotomy, and
the model's **structural choice ②** (two separable injury channels) comes from here.

---

## 2. Haematoma expansion — the model's central integrated quantity

14. Brott T, Broderick J, Kothari R, et al. **Early hemorrhage growth in patients with
    intracerebral hemorrhage.** *Stroke.* 1997;28:1-5.
    <https://pubmed.ncbi.nlm.nih.gov/8996478/>
15. Davis SM, Broderick J, Hennerici M, et al. **Hematoma growth is a determinant of
    mortality and poor outcome after intracerebral hemorrhage.** *Neurology.*
    2006;66:1175-1181. <https://pubmed.ncbi.nlm.nih.gov/16636233/>
16. Al-Shahi Salman R, Frantzias J, Lee RJ, et al. **Absolute risk and predictors of the
    growth of acute spontaneous intracerebral haemorrhage: a systematic review and
    meta-analysis of individual patient data.** *Lancet Neurol.* 2018;17:885-894.
    <https://pubmed.ncbi.nlm.nih.gov/30120039/>
17. Morotti A, Boulouis G, Dowlatshahi D, et al. **Intracerebral haemorrhage expansion:
    definitions, predictors, and prevention.** *Lancet Neurol.* 2023;22:159-171.
    <https://pubmed.ncbi.nlm.nih.gov/36309041/>
18. Kothari RU, Brott T, Broderick JP, et al. **The ABCs of measuring intracerebral
    hemorrhage volumes.** *Stroke.* 1996;27:1304-1305.
    <https://pubmed.ncbi.nlm.nih.gov/8711791/>
19. Wada R, Aviv RI, Fox AJ, et al. **CT angiography "spot sign" predicts hematoma
    expansion in acute intracerebral hemorrhage.** *Stroke.* 2007;38:1257-1262.
    <https://pubmed.ncbi.nlm.nih.gov/17322083/>
20. Demchuk AM, Dowlatshahi D, Rodriguez-Luna D, et al. **Prediction of haematoma growth
    and outcome in patients with intracerebral haemorrhage using the CT-angiography spot
    sign (PREDICT): a prospective observational study.** *Lancet Neurol.* 2012;11:307-314.
    <https://pubmed.ncbi.nlm.nih.gov/22405630/>
21. Delgado Almandoz JE, Yoo AJ, Stone MJ, et al. **The spot sign score in primary
    intracerebral hemorrhage identifies patients at highest risk of in-hospital mortality
    and poor outcome among survivors.** *Stroke.* 2010;41:54-60.
    <https://pubmed.ncbi.nlm.nih.gov/19910545/>
22. Li Q, Zhang G, Huang YJ, et al. **Blend Sign on Computed Tomography: Novel and
    Reliable Predictor for Early Hematoma Growth in Patients With Intracerebral
    Hemorrhage.** *Stroke.* 2015;46:2119-2123.
    <https://pubmed.ncbi.nlm.nih.gov/26089330/>
23. Li Q, Zhang G, Xiong X, et al. **Black Hole Sign: Novel Imaging Marker That Predicts
    Hematoma Growth in Patients With Intracerebral Hemorrhage.** *Stroke.* 2016;47:1777-1781.
    <https://pubmed.ncbi.nlm.nih.gov/27174523/>
24. Blacquiere D, Demchuk AM, Al-Hazzaa M, et al. **Intracerebral Hematoma Morphologic
    Appearance on Noncontrast Computed Tomography Predicts Significant Hematoma
    Expansion.** *Stroke.* 2015;46:3111-3116.
    <https://pubmed.ncbi.nlm.nih.gov/26451019/>

**▶ Model link** — this group of references is the quantitative target of **structural choice ①**. Refs 14-16
supply the untreated expansion rate of 20-38% and the fact that expansion mostly occurs in the first few hours, and that is
the basis for calibrating `KBLEED` (0.085) and `KSEAL` (0.10, which stretches the bleeding window to 6-12 hours).
In the draft, when `KSEAL = 0.80`, bleeding terminated within 1.5 hours and so was over before the antihypertensive arrived,
and the BP effect came out far smaller than in the literature — the model was fixed because of
this group of references. Ref 15 is the basis of the `dV24 -> NIHSS -> mRS` path, and refs 19-24 correspond to the map's
CTA spot sign / blend sign / black hole sign nodes and the `SpotSign` state.

---

## 3. Blood pressure management — two pathways with opposite signs

25. Anderson CS, Heeley E, Huang Y, et al. **Rapid blood-pressure lowering in patients
    with acute intracerebral hemorrhage (INTERACT2).** *N Engl J Med.* 2013;368:2355-2365.
    <https://pubmed.ncbi.nlm.nih.gov/23713578/>
26. Qureshi AI, Palesch YY, Barsan WG, et al. **Intensive Blood-Pressure Lowering in
    Patients with Acute Cerebral Hemorrhage (ATACH-2).** *N Engl J Med.* 2016;375:1033-1043.
    <https://pubmed.ncbi.nlm.nih.gov/27276234/>
27. Ma L, Hu X, Song L, et al. **The third Intensive Care Bundle with Blood Pressure
    Reduction in Acute Cerebral Haemorrhage Trial (INTERACT3): an international,
    stepped wedge cluster randomised controlled trial.** *Lancet.* 2023;402:27-40.
    <https://pubmed.ncbi.nlm.nih.gov/37245517/>
28. Li G, Lin Y, Yang J, et al. **Intensive Ambulance-Delivered Blood-Pressure Reduction
    in Hyperacute Stroke (INTERACT4).** *N Engl J Med.* 2024;390:1862-1872.
    <https://pubmed.ncbi.nlm.nih.gov/38752650/>
29. Anderson CS, Huang Y, Wang JG, et al. **Intensive blood pressure reduction in acute
    cerebral haemorrhage trial (INTERACT): a randomised pilot trial.** *Lancet Neurol.*
    2008;7:391-399. <https://pubmed.ncbi.nlm.nih.gov/18396107/>
30. Moullaali TJ, Wang X, Martin RH, et al. **Blood pressure control and clinical outcomes
    in acute intracerebral haemorrhage: a preplanned pooled analysis of individual
    participant data.** *Lancet Neurol.* 2019;18:857-864.
    <https://pubmed.ncbi.nlm.nih.gov/31397290/>
31. Manning L, Hirakawa Y, Arima H, et al. **Blood pressure variability and outcome after
    acute intracerebral haemorrhage: a post-hoc analysis of INTERACT2.** *Lancet Neurol.*
    2014;13:364-373. <https://pubmed.ncbi.nlm.nih.gov/24530176/>
32. Qureshi AI, Huang W, Lobanova I, et al. **Outcomes of Intensive Systolic Blood
    Pressure Reduction in Patients With Intracerebral Hemorrhage and Excessively High
    Initial Systolic Blood Pressure: Post Hoc Analysis of a Randomized Clinical Trial.**
    *JAMA Neurol.* 2020;77:1355-1365. <https://pubmed.ncbi.nlm.nih.gov/32897310/>
33. Butcher K, Jeerakathil T, Emery D, et al. **The Intracerebral Haemorrhage Acutely
    Decreasing Arterial Pressure Trial: ICH ADAPT.** *Int J Stroke.* 2010;5:227-233.
    <https://pubmed.ncbi.nlm.nih.gov/20536619/>
34. Gould B, McCourt R, Asdaghi N, et al. **Autoregulation of cerebral blood flow is
    preserved in primary intracerebral hemorrhage.** *Stroke.* 2013;44:1726-1728.
    <https://pubmed.ncbi.nlm.nih.gov/23619129/>
35. Toyoda K, Koga M, Yamamoto H, et al. **Intensive blood pressure lowering with
    nicardipine and outcomes after intracerebral hemorrhage: An individual participant
    data systematic review.** *Int J Stroke.* 2022;17:494-505.
    <https://pubmed.ncbi.nlm.nih.gov/34542358/>

**▶ Model link** — both the basis for and the refutation of **structural choice ③**. The effect sizes supplied by refs 25-26
(INTERACT2 mean growth 4.5 vs 5.3 mL, ATACH-2 expansion 18.9% vs 24.4%) set the upper bound on
the calibration of `KBLEED` · `EMXNIC` · `EC5NIC`. In the draft the model had been tuned to produce a 40% reduction,
but that is **reproducing a result the trials did not obtain**, so the target was moved back towards
the literature. Ref 32 is the basis of `TCPP` (cumulative time with CPP<60) · the renal term and of scenario 4 (excessive
BP lowering = minimum bleeding + worst outcome), and refs 33-34 are the basis of `AUTOR` · `CBFREG` · `RSHIFT` (rightward
shift of hypertensive autoregulation). Ref 31 corresponds to the `SBPV` node. Ref 27 supports the multiplicative structure of
scenario 6 (care bundle), and ref 28 the map's `Prehospital` node.

---

## 4. Haemostatic therapy · anticoagulation reversal

36. Mayer SA, Brun NC, Begtrup K, et al. **Efficacy and safety of recombinant activated
    factor VII for acute intracerebral hemorrhage (FAST).** *N Engl J Med.*
    2008;358:2127-2137. <https://pubmed.ncbi.nlm.nih.gov/18480205/>
37. Sprigg N, Flaherty K, Appleton JP, et al. **Tranexamic acid for hyperacute primary
    IntraCerebral Haemorrhage (TICH-2): an international randomised, placebo-controlled,
    phase 3 superiority trial.** *Lancet.* 2018;391:2107-2115.
    <https://pubmed.ncbi.nlm.nih.gov/29778325/>
38. Baharoglu MI, Cordonnier C, Al-Shahi Salman R, et al. **Platelet transfusion versus
    standard care after acute stroke due to spontaneous cerebral haemorrhage associated
    with antiplatelet therapy (PATCH): a randomised, open-label, phase 3 trial.**
    *Lancet.* 2016;387:2605-2613. <https://pubmed.ncbi.nlm.nih.gov/27178479/>
39. Steiner T, Poli S, Griebe M, et al. **Fresh frozen plasma versus prothrombin complex
    concentrate in patients with intracranial haemorrhage related to vitamin K
    antagonists (INCH): a randomised trial.** *Lancet Neurol.* 2016;15:566-573.
    <https://pubmed.ncbi.nlm.nih.gov/27302126/>
40. Connolly SJ, Sharma M, Cohen AT, et al. **Andexanet for Factor Xa Inhibitor-Associated
    Acute Intracerebral Hemorrhage (ANNEXA-I).** *N Engl J Med.* 2024;390:1745-1755.
    <https://pubmed.ncbi.nlm.nih.gov/38749032/>
41. Connolly SJ, Crowther M, Eikelboom JW, et al. **Full Study Report of Andexanet Alfa
    for Bleeding Associated with Factor Xa Inhibitors (ANNEXA-4).** *N Engl J Med.*
    2019;380:1326-1335. <https://pubmed.ncbi.nlm.nih.gov/30730782/>
42. Pollack CV Jr, Reilly PA, van Ryn J, et al. **Idarucizumab for Dabigatran Reversal —
    Full Cohort Analysis (RE-VERSE AD).** *N Engl J Med.* 2017;377:431-441.
    <https://pubmed.ncbi.nlm.nih.gov/28693366/>
43. Siegal DM, Curnutte JT, Connolly SJ, et al. **Andexanet Alfa for the Reversal of
    Factor Xa Inhibitor Activity.** *N Engl J Med.* 2015;373:2413-2424.
    <https://pubmed.ncbi.nlm.nih.gov/26559317/>
44. Flaherty ML, Tao H, Haverbusch M, et al. **Warfarin use leads to larger intracerebral
    hematomas.** *Neurology.* 2008;71:1084-1089.
    <https://pubmed.ncbi.nlm.nih.gov/18824672/>
45. Kuramatsu JB, Gerner ST, Schellinger PD, et al. **Anticoagulant reversal, blood
    pressure levels, and anticoagulant resumption in patients with anticoagulation-related
    intracerebral hemorrhage.** *JAMA.* 2015;313:824-836.
    <https://pubmed.ncbi.nlm.nih.gov/25710659/>
46. Eddleston M, de la Torre JC, Oldstone MB, Loskutoff DJ, Edgington TS, Mackman N.
    **Astrocytes are the primary source of tissue factor in the murine central nervous
    system.** *J Clin Invest.* 1993;92:349-358.
    <https://pubmed.ncbi.nlm.nih.gov/8326003/>

**▶ Model link** — the whole of clusters 3-5 and the `CLOT` · `FII` · `PLAS` · `PLTF` equations.
Ref 46 is the evidence that the brain expresses tissue factor (TF) very abundantly, and hence the reason for the `TF_brain -> FVIIa_TF`
path and for `KFORM` taking a large value (3.0). Refs 44-45 support the `FII0 = 20%` (INR 3.0)
scenario and the longer, larger expansion of OAC-ICH, but **quantitatively it is about twofold** —
in the draft, using a single multiplicative term (`THRGEN*PLTF`) gave 3.5-fold, and so the `CLOT` equation was
split into a platelet arm (`WPLT`) and a coagulation arm (`WFIB`). That split is also consistent with ref 38 (PATCH: platelet
transfusion is **harmful** in antiplatelet-associated ICH), and it is the reason platelet transfusion and DDAVP have
no effect at all on VKA bleeding in the model. Refs 40 and 43 are the basis of `KBINDA` · `KOFFA` · `KELCPA` and of the
**rebound** after andexanet is stopped (which emerges in the model from dissociation + redistribution, with no separate term),
and ref 36, where rFVIIa reduced expansion but did not improve outcome (and increased thrombosis), is the
counter-example showing **the limits of a single-factor intervention**.

---

## 5. Intracranial pressure · perfusion · autoregulation

47. Mokri B. **The Monro-Kellie hypothesis: applications in CSF volume depletion.**
    *Neurology.* 2001;56:1746-1748. <https://pubmed.ncbi.nlm.nih.gov/11425944/>
48. Marmarou A, Shulman K, LaMorgese J. **Compartmental analysis of compliance and outflow
    resistance of the cerebrospinal fluid system.** *J Neurosurg.* 1975;43:523-534.
    <https://pubmed.ncbi.nlm.nih.gov/1181384/>
49. Davson H, Hollingsworth G, Segal MB. **The mechanism of drainage of the cerebrospinal
    fluid.** *Brain.* 1970;93:665-678. <https://pubmed.ncbi.nlm.nih.gov/5490270/>
50. Zazulia AR, Diringer MN, Videen TO, et al. **Hypoperfusion without ischemia
    surrounding acute intracerebral hemorrhage.** *J Cereb Blood Flow Metab.*
    2001;21:804-810. <https://pubmed.ncbi.nlm.nih.gov/11435792/>
51. Oeinck M, Neunhoeffer F, Buttler KJ, et al. **Dynamic cerebral autoregulation in acute
    intracerebral hemorrhage.** *Stroke.* 2013;44:2722-2728.
    <https://pubmed.ncbi.nlm.nih.gov/23943213/>
52. Reinhard M, Neunhoeffer F, Gerds TA, et al. **Secondary decline of cerebral
    autoregulation is associated with worse outcome after intracerebral hemorrhage.**
    *Intensive Care Med.* 2010;36:264-271.
    <https://pubmed.ncbi.nlm.nih.gov/19838669/>
53. Czosnyka M, Smielewski P, Kirkpatrick P, Laing RJ, Menon D, Pickard JD. **Continuous
    assessment of the cerebral vasomotor reactivity in head injury.** *Neurosurgery.*
    1997;41:11-17. <https://pubmed.ncbi.nlm.nih.gov/9218290/>

**▶ Model link** — ref 47 is the basis of the design that treats the `MonroKellie` identity and `VADD` as **the residual of the
volume balance rather than as a state**. Refs 48-49 (the Davson equation: absorption = (ICP − Pss) × conductance) are
the source of `KCSFC = 4.2` and `PSS = 5`, and these two values balance production (21 mL/h) exactly in the
undamaged state — in the draft this balance was absent and **the CSF diverged negative**. `VCSFMIN`
(18 mL of expressible CSF) makes the pressure-volume curve **biphasic** rather than exponential,
which is the real physiology of exhausted compensation. Ref 50 is the basis of `CBFISC = 0.70` —
perihaematomal tissue is already oligaemic, so the threshold at which injury begins has to lie **higher** than the classical
infarction threshold, and when `CBFISC = 0.55` the perfusion channel did not operate at all,
so the U-shaped curve failed to emerge. Refs 51-53 are the basis of the dynamic decline of `AUTOR`, the `PRx` node,
and the `KISCH` term by which loss of autoregulation worsens the outcome.

---

## 6. Perihaematomal oedema

54. Wagner KR, Xi G, Hua Y, et al. **Lobar intracerebral hemorrhage model in pigs: rapid
    edema development in perihematomal white matter.** *Stroke.* 1996;27:490-497.
    <https://pubmed.ncbi.nlm.nih.gov/8610319/>
55. Xi G, Wagner KR, Keep RF, et al. **Role of blood clot formation on early edema
    development after experimental intracerebral hemorrhage.** *Stroke.* 1998;29:2580-2586.
    <https://pubmed.ncbi.nlm.nih.gov/9836771/>
56. Hua Y, Keep RF, Hoff JT, Xi G. **Thrombin preconditioning attenuates brain edema
    induced by erythrocytes and iron.** *J Cereb Blood Flow Metab.* 2003;23:1448-1454.
    <https://pubmed.ncbi.nlm.nih.gov/14663340/>
57. Urday S, Kimberly WT, Beslow LA, et al. **Targeting secondary injury in intracerebral
    haemorrhage — perihaematomal oedema.** *Nat Rev Neurol.* 2015;11:111-122.
    <https://pubmed.ncbi.nlm.nih.gov/25623787/>
58. Volbers B, Giede-Jeppe A, Gerner ST, et al. **Peak perihemorrhagic edema correlates
    with functional outcome in intracerebral hemorrhage.** *Neurology.* 2018;90:e1005-e1012.
    <https://pubmed.ncbi.nlm.nih.gov/29453243/>
59. Simard JM, Chen M, Tarasov KV, Bhatta S, Ivanova S, Melnitchenko L, Tsymbalyuk N,
    West GA, Gerzanich V. **Newly expressed SUR1-regulated NC(Ca-ATP) channel mediates
    cerebral edema after ischemic stroke.** *Nat Med.* 2006;12:433-440.
    <https://pubmed.ncbi.nlm.nih.gov/16550187/>
60. Jiang B, Li L, Chen Q, et al. **Glibenclamide Attenuates Neuroinflammation and
    Promotes Neurological Recovery After Intracerebral Hemorrhage in Aged Rats.**
    *Front Aging Neurosci.* 2021;13:729652.
    <https://pubmed.ncbi.nlm.nih.gov/34512312/>
61. Chen Y, Chen S, Chang J, Wei J, Feng M, Wang R. **Perihematomal Edema After
    Intracerebral Hemorrhage: An Update on Pathogenesis, Risk Factors, and Therapeutic
    Advances.** *Front Immunol.* 2021;12:740632.
    <https://pubmed.ncbi.nlm.nih.gov/34737745/>
62. Gong Y, Xi G, Hu H, et al. **Complement inhibition attenuates brain edema and
    neurological deficits induced by thrombin.** *Acta Neurochir Suppl.* 2005;95:389-392.
    <https://pubmed.ncbi.nlm.nih.gov/16463887/>

**▶ Model link** — the **separation into two compartments**, `OEDE` (ionic · vasogenic, peaking at 3-5 days) and `OEDL`
(cytotoxic · iron-driven, 7-14 days), comes from refs 54-57. Refs 55-56 and 62 are the basis of the `THRT -> PAR1 -> BBBP -> OEDE`
path and of `KTHR` · `KTHRI`, and refs 59-60 are the `SUR1_TRPM4` node and the `KSUR1` term.
Ref 58 is the evidence that the oedema peak predicts outcome **independently** of haematoma volume, and hence the reason the `PHE_rel`
node and `STRAIN` include `OEDE + OEDL` as well as `VHEM` alone.
In model validation the total PHE peak is 26.2 mL (4.3 days), agreeing with the 3-5 day peak in the literature.

---

## 7. Haem · iron · ferroptosis — the delayed chemical channel

63. Wu J, Hua Y, Keep RF, Nakamura T, Hoff JT, Xi G. **Iron and iron-handling proteins in
    the brain after intracerebral hemorrhage.** *Stroke.* 2003;34:2964-2969.
    <https://pubmed.ncbi.nlm.nih.gov/14615611/>
64. Pérez de la Ossa N, Sobrino T, Silva Y, et al. **Iron-related brain damage in patients
    with intracerebral hemorrhage.** *Stroke.* 2010;41:810-813.
    <https://pubmed.ncbi.nlm.nih.gov/20185788/>
65. Hua Y, Keep RF, Hoff JT, Xi G. **Deferoxamine therapy for intracerebral hemorrhage.**
    *Acta Neurochir Suppl.* 2008;105:3-6.
    <https://pubmed.ncbi.nlm.nih.gov/19066072/>
66. Hu S, Hua Y, Keep RF, Feng H, Xi G. **Deferoxamine therapy reduces brain hemin
    accumulation after intracerebral hemorrhage in piglets.** *Exp Neurol.* 2019;318:244-250.
    <https://pubmed.ncbi.nlm.nih.gov/31078524/>
67. Selim M, Foster LD, Moy CS, et al. **Deferoxamine mesylate in patients with
    intracerebral haemorrhage (i-DEF): a multicentre, randomised, placebo-controlled,
    double-blind phase 2 trial.** *Lancet Neurol.* 2019;18:428-438.
    <https://pubmed.ncbi.nlm.nih.gov/30898550/>
68. Selim M, Yeatts S, Goldstein JN, et al. **Safety and tolerability of deferoxamine
    mesylate in patients with acute intracerebral hemorrhage.** *Stroke.* 2011;42:3067-3074.
    <https://pubmed.ncbi.nlm.nih.gov/21868742/>
69. Zeng L, Tan L, Li H, Zhang Q, Li Y, Guo J. **Deferoxamine therapy for intracerebral
    hemorrhage: A systematic review.** *PLoS One.* 2018;13:e0193615.
    <https://pubmed.ncbi.nlm.nih.gov/29566000/>
70. Zille M, Karuppagounder SS, Chen Y, et al. **Neuronal Death After Hemorrhagic Stroke
    In Vitro and In Vivo Shares Features of Ferroptosis and Necroptosis.** *Stroke.*
    2017;48:1033-1043. <https://pubmed.ncbi.nlm.nih.gov/28250197/>
71. Dixon SJ, Lemberg KM, Lamprecht MR, et al. **Ferroptosis: an iron-dependent form of
    nonapoptotic cell death.** *Cell.* 2012;149:1060-1072.
    <https://pubmed.ncbi.nlm.nih.gov/22632970/>
72. Zhao X, Song S, Sun G, et al. **Neuroprotective role of haptoglobin after
    intracerebral hemorrhage.** *J Neurosci.* 2009;29:15819-15827.
    <https://pubmed.ncbi.nlm.nih.gov/20016097/>
73. Chen-Roetling J, Regan RF. **Hemopexin increases the neurotoxicity of hemoglobin when
    haptoglobin is absent.** *J Neurochem.* 2018;145:464-473.
    <https://pubmed.ncbi.nlm.nih.gov/29500821/>
74. Wang Y, Kinzie E, Berger FG, Lim SK, Baumann H. **Haptoglobin, an
    inflammation-inducible plasma protein.** *Redox Rep.* 2001;6:379-385.
    <https://pubmed.ncbi.nlm.nih.gov/11865981/>
75. Zhao X, Sun G, Zhang J, et al. **Transcription factor Nrf2 protects the brain from
    damage produced by intracerebral hemorrhage.** *Stroke.* 2007;38:3280-3286.
    <https://pubmed.ncbi.nlm.nih.gov/17962605/>
76. Summers MR, Jacobs A, Tudway D, Perera P, Ricketts C. **Studies in desferrioxamine and
    ferrioxamine metabolism in normal and iron-loaded subjects.** *Br J Haematol.*
    1979;42:547-555. <https://pubmed.ncbi.nlm.nih.gov/476006/>

**▶ Model link** — the whole chemical channel of **structural choice ②**:
`RBC_lysis -> HEME -> HO1 -> FEII -> Fenton -> LPO -> ferroptosis -> NEUR`.
Refs 63 and 64 are the basis of `KHEMI` · `YFE` · `FERR` (ferritin sequestration) and the serum ferritin node,
and refs 70-71 are the mechanistic basis of the `KFERRO2` (lipid peroxidation → neuronal death) term.
Refs 72-74 are the `HPG` parameter (scavenging capacity by haptoglobin genotype) and the `KSCAV` term,
and ref 75 is the `Nrf2 -> HO1/FERR/MG2` path.
Refs 65-69 and 76 are the basis of deferoxamine PK (`CLDFO`, t½ ≈ 18 min) and of `KCHEL`, and
**in the model DFO has not a single edge going to `VHEM`** — on validation, the i-DEF scenario came out with a
24-hour haematoma volume **exactly identical** to the control while only the 14-day iron exposure AUC fell 68.8 → 48.1
(−30%), and setting `KCHEL = 0` makes that effect disappear entirely. This is how the model structurally reproduces
the "mRS shift without volume change" pattern that ref 67 reported.

---

## 8. Neuroinflammation · haematoma clearance

77. Mracsko E, Veltkamp R. **Neuroinflammation after intracerebral hemorrhage.**
    *Front Cell Neurosci.* 2014;8:388. <https://pubmed.ncbi.nlm.nih.gov/25477782/>
78. Ma Q, Chen S, Hu Q, Feng H, Zhang JH, Tang J. **NLRP3 inflammasome contributes to
    inflammation after intracerebral hemorrhage.** *Ann Neurol.* 2014;75:209-219.
    <https://pubmed.ncbi.nlm.nih.gov/24273204/>
79. Zhao X, Sun G, Zhang J, et al. **Hematoma resolution as a target for intracerebral
    hemorrhage treatment: role for peroxisome proliferator-activated receptor gamma in
    microglia/macrophages.** *Ann Neurol.* 2007;61:352-362.
    <https://pubmed.ncbi.nlm.nih.gov/17457822/>
80. Gonzales NR, Shah J, Sangha N, et al. **Design of a prospective, dose-escalation study
    evaluating the Safety of Pioglitazone for Hematoma Resolution in Intracerebral
    Hemorrhage (SHRINC).** *Int J Stroke.* 2013;8:388-396.
    <https://pubmed.ncbi.nlm.nih.gov/22340518/>
81. Chang CF, Massey J, Osherov A, Angenendt da Costa LH, Sansing LH. **Bexarotene
    Enhances Macrophage Erythrophagocytosis and Hematoma Clearance in Experimental
    Intracerebral Hemorrhage.** *Stroke.* 2020;51:612-618.
    <https://pubmed.ncbi.nlm.nih.gov/31826730/>
82. Fu Y, Hao J, Zhang N, et al. **Fingolimod for the treatment of intracerebral
    hemorrhage: a 2-arm proof-of-concept study.** *JAMA Neurol.* 2014;71:1092-1101.
    <https://pubmed.ncbi.nlm.nih.gov/25003359/>

**▶ Model link** — the `MG1` (pro-inflammatory) · `MG2` (reparative · phagocytic) · `NEU` · `IL6` · `IL10` equations.
Refs 79 and 81 are the basis for making the haematoma resorption rate in the model not a constant but **`KRESE = KRES0 × MG2`**,
that is, rate-limited by reparative microglia (so an anti-inflammatory intervention may also
slow clearance). Ref 78 is the `NLRP3 -> SUR1_TRPM4` · `IL6` path,
and ref 82 is the `Fingolimod` node.

**One modelling caveat** — during validation, writing the `IL6`↔`NEU` feedback linearly
took the loop gain above 1 and the cytokine states **diverged** (IL-6 running away to 18-fold at 24 hours).
Only when both directions are changed to saturating (Michaelis) functions does the inflammatory response become the **finite,
self-resolving pulse** the literature describes. This is not a cosmetic fix but a structural requirement.

---

## 9. Intraventricular haemorrhage · hydrocephalus

83. Hanley DF, Lane K, McBee N, et al. **Thrombolytic removal of intraventricular
    haemorrhage in treatment of severe stroke: results of the randomised, multicentre,
    multiregion, placebo-controlled CLEAR III trial.** *Lancet.* 2017;389:603-611.
    <https://pubmed.ncbi.nlm.nih.gov/28081952/>
84. Graeb DA, Robertson WD, Lapointe JS, Nugent RA, Harrison PB. **Computed tomographic
    diagnosis of intraventricular hemorrhage. Etiology and prognosis.** *Radiology.*
    1982;143:91-96. <https://pubmed.ncbi.nlm.nih.gov/6977795/>
85. Hallevi H, Albright KC, Aronowski J, et al. **Intraventricular hemorrhage: Anatomic
    relationships and clinical implications.** *Neurology.* 2008;70:848-852.
    <https://pubmed.ncbi.nlm.nih.gov/18332342/>
86. Roh DJ, Boehme A, Mamoon R, et al. **Intraventricular Hemorrhage Expansion in the
    CLEAR III Trial: A Post Hoc Exploratory Analysis.** *Stroke.* 2022;53:1847-1853.
    <https://pubmed.ncbi.nlm.nih.gov/35086362/>
87. Strahle JM, Garton HJ, Maher CO, Muraszko KM, Keep RF, Xi G. **Role of hemoglobin and
    iron in hydrocephalus after neonatal intraventricular hemorrhage.** *Neurosurgery.*
    2014;75:696-705. <https://pubmed.ncbi.nlm.nih.gov/25121790/>

**▶ Model link** — the `VIVH` · `VCSF` · `EVD` · `IVTPA` equations. Ref 83 is the basis of `EMXTPA` · `KIVHCL` and of
scenario 12, and **the central asymmetry of CLEAR III (mortality falls but mRS 0-3 is almost
unchanged)** is reproduced in this model too. Refs 84-85 are the `IVH_cast` · Graeb nodes,
and ref 87 is the `Ependymal` · `ShuntDep` path by which iron after IVH contributes to hydrocephalus.
The **asymmetry** that the blockade of CSF absorption (`1 - VIVH/(KIVHB+VIVH)`) is not applied to production is
the whole of the IVH → hydrocephalus mechanism.

---

## 10. Surgical strategies

88. Mendelow AD, Gregson BA, Fernandes HM, et al. **Early surgery versus initial
    conservative treatment in patients with spontaneous supratentorial intracerebral
    haematomas in the International Surgical Trial in Intracerebral Haemorrhage (STICH):
    a randomised trial.** *Lancet.* 2005;365:387-397.
    <https://pubmed.ncbi.nlm.nih.gov/15680453/>
89. Mendelow AD, Gregson BA, Rowan EN, et al. **Early surgery versus initial conservative
    treatment in patients with spontaneous supratentorial lobar intracerebral haematomas
    (STICH II): a randomised trial.** *Lancet.* 2013;382:397-408.
    <https://pubmed.ncbi.nlm.nih.gov/23726393/>
90. Hanley DF, Thompson RE, Rosenblum M, et al. **Efficacy and safety of minimally
    invasive surgery with thrombolysis in intracerebral haemorrhage evacuation
    (MISTIE III): a randomised, controlled, open-label, blinded endpoint phase 3 trial.**
    *Lancet.* 2019;393:1021-1032. <https://pubmed.ncbi.nlm.nih.gov/30739747/>
91. Pradilla G, Ratcliff JJ, Hall AJ, et al. **Trial of Early Minimally Invasive Removal
    of Intracerebral Hemorrhage (ENRICH).** *N Engl J Med.* 2024;390:1277-1289.
    <https://pubmed.ncbi.nlm.nih.gov/38598795/>
92. Beck J, Fung C, Strbian D, et al. **Decompressive craniectomy plus best medical
    treatment versus best medical treatment alone for spontaneous severe deep
    supratentorial intracerebral haemorrhage: a randomised controlled clinical trial**
    (SWITCH). *Lancet.* 2024;403:2395-2404.
    <https://pubmed.ncbi.nlm.nih.gov/38761811/>

**▶ Model link** — ref 90 is one of the most important structural lessons of the model: **the mediator of the
effect is not the surgical approach but `EOT_vol` (the residual volume at the end of treatment)**, and the benefit in MISTIE III
was confined to the subgroup achieving ≤15 mL. In the model, therefore, surgery reduces `VHEM` through `FEVAC`
while at the same time raising `NOPEN` through `SURGREB` (procedure-related rebleeding);
it acts on the mass-effect channel only and cannot touch iron already liberated. Refs 88-89 (the neutral results of
STICH I/II) are the counter-example that "the volume was reduced" alone does not improve outcome, and hence the reason the `STRTHR`
term exists (below-threshold strain makes mass injury negligible). Ref 92
corresponds to `DecompCranx -> VDC` (the compensatory volume gained when the skull is opened).

---

## 11. Secondary prevention · cerebral amyloid angiopathy

93. SPS3 Study Group. **Blood-pressure targets in patients with recent lacunar stroke:
    the SPS3 randomised trial.** *Lancet.* 2013;382:507-515.
    <https://pubmed.ncbi.nlm.nih.gov/23726159/>
94. PROGRESS Collaborative Group. **Randomised trial of a perindopril-based
    blood-pressure-lowering regimen among 6105 individuals with previous stroke or
    transient ischaemic attack.** *Lancet.* 2001;358:1033-1041.
    <https://pubmed.ncbi.nlm.nih.gov/11589932/>
95. RESTART Collaboration. **Effects of antiplatelet therapy after stroke due to
    intracerebral haemorrhage (RESTART): a randomised, open-label trial.** *Lancet.*
    2019;393:2613-2623. <https://pubmed.ncbi.nlm.nih.gov/31128924/>
96. SoSTART Collaboration. **Effects of oral anticoagulation for atrial fibrillation after
    spontaneous intracranial haemorrhage in the UK: a randomised, open-label, assessor-
    masked, pilot-phase, non-inferiority trial.** *Lancet Neurol.* 2021;20:842-853.
    <https://pubmed.ncbi.nlm.nih.gov/34487722/>
97. Schreuder FHBM, van Nieuwenhuizen KM, Hofmeijer J, et al. **Apixaban versus no
    anticoagulation after anticoagulation-associated intracerebral haemorrhage in patients
    with atrial fibrillation (APACHE-AF): a randomised, open-label, phase 2 trial.**
    *Lancet Neurol.* 2021;20:907-916. <https://pubmed.ncbi.nlm.nih.gov/34687635/>
98. Charidimou A, Boulouis G, Frosch MP, et al. **The Boston criteria version 2.0 for
    cerebral amyloid angiopathy: a multicentre, retrospective, MRI-neuropathology
    diagnostic accuracy study.** *Lancet Neurol.* 2022;21:714-725.
    <https://pubmed.ncbi.nlm.nih.gov/35841910/>
99. Charidimou A, Boulouis G, Gurol ME, et al. **Emerging concepts in sporadic cerebral
    amyloid angiopathy.** *Brain.* 2017;140:1829-1850.
    <https://pubmed.ncbi.nlm.nih.gov/28334869/>
100. Charidimou A, Krishnan A, Werring DJ, Rolf Jäger H. **Cerebral microbleeds: a guide to
     detection and clinical relevance in different disease settings.** *Neuroradiology.*
     2013;55:655-674. <https://pubmed.ncbi.nlm.nih.gov/23708941/>

**▶ Model link** — the whole of cluster 17 (secondary prevention) in the map, and the `Recurrence` node's
risk distinction between deep (about 2%/year) and lobar · CAA (about 7%/year). Refs 98-100 are the `CAA` · `Microbleeds` ·
`cSS` · `APOE` nodes and the Boston v2.0 criteria, refs 93-94 are `BP_longterm`,
and refs 95-97 are the `OAC_restart` · `Antiplatelet_restart` dilemma.
(This cluster is included in the map but, being outside the 90-day simulation horizon, is not
implemented in the ODEs — recurrence is not a state variable of the model.)

---

## 12. Prognosis · biomarkers

101. Hemphill JC 3rd, Bonovich DC, Besmertis L, Manley GT, Johnston SC. **The ICH score: a
     simple, reliable grading scale for intracerebral hemorrhage.** *Stroke.*
     2001;32:891-897. <https://pubmed.ncbi.nlm.nih.gov/11283388/>
102. Sembill JA, Castello JP, Sprügel MI, et al. **Multicenter Validation of the max-ICH
     Score in Intracerebral Hemorrhage.** *Ann Neurol.* 2021;89:474-484.
     <https://pubmed.ncbi.nlm.nih.gov/33222266/>
103. Rost NS, Smith EE, Chang Y, et al. **Prediction of functional outcome in patients with
     primary intracerebral hemorrhage: the FUNC score.** *Stroke.* 2008;39:2304-2309.
     <https://pubmed.ncbi.nlm.nih.gov/18556582/>
104. Becker KJ, Baxter AB, Cohen WA, et al. **Withdrawal of support in intracerebral
     hemorrhage may lead to self-fulfilling prophecies.** *Neurology.* 2001;56:766-772.
     <https://pubmed.ncbi.nlm.nih.gov/11274312/>
105. Foerch C, Niessner M, Back T, et al. **Diagnostic accuracy of plasma glial fibrillary
     acidic protein for differentiating intracerebral hemorrhage and cerebral ischemia in
     patients with symptoms of acute stroke.** *Clin Chem.* 2012;58:237-245.
     <https://pubmed.ncbi.nlm.nih.gov/22125303/>
106. Wang X, Moullaali TJ, Li Q, et al. **Utility-Weighted Modified Rankin Scale Scores for
     the Assessment of Stroke Outcome: Pooled Analysis of 20 000+ Patients.** *Stroke.*
     2020;51:2411-2417. <https://pubmed.ncbi.nlm.nih.gov/32640944/>

**▶ Model link** — `ICHSC` in `$TABLE` (exactly the original scoring rule of ref 101: GCS · volume ≥30 mL ·
presence of IVH · infratentorial origin · age ≥80), the `PMORT` logit (`MORTA` · `MORTB`), `UWMRS` (the utility
weights of ref 106), and the `GFAP_ratio` node (ref 105). Ref 104 is included explicitly as the map's `DNR_bias`
node — **a self-fulfilling prophecy that confounds the interpretation of the results of every ICH trial**, and one that
must be allowed for when comparing the "treatment effect" the model reproduces with actual trial results.

**A caution on calibration** — `ICHSC` is by design a **baseline score**. In the draft the ICH score was computed from
the state at 90 days to derive mortality, which is scoring on the basis of a patient who has recovered and so
inverts the intent. It is now combined from the 24-hour score + the cumulative intracranial pressure exposure
(`TICP`).

---

## 13. Pharmacokinetics · pharmacodynamics · QSP methodology

107. Cook E, Clifton GG, Vargas R, et al. **Pharmacokinetics, pharmacodynamics, and
     minimum effective clinical dose of intravenous nicardipine.** *Clin Pharmacol Ther.*
     1990;47:706-718. <https://pubmed.ncbi.nlm.nih.gov/2357865/>
108. Keating GM. **Clevidipine: a review of its use for managing blood pressure in
     perioperative and intensive care settings.** *Drugs.* 2014;74:1947-1960.
     <https://pubmed.ncbi.nlm.nih.gov/25312594/>
109. Li S, Ahmadzia HK, Guo D, et al. **Population pharmacokinetics and pharmacodynamics
     of tranexamic acid in women undergoing caesarean delivery.** *Br J Clin Pharmacol.*
     2021;87:3579-3589. <https://pubmed.ncbi.nlm.nih.gov/33576009/>
110. Elmokadem A, Riggs MM, Baron KT. **Quantitative Systems Pharmacology and
     Physiologically-Based Pharmacokinetic Modeling With mrgsolve: A Hands-On Tutorial.**
     *CPT Pharmacometrics Syst Pharmacol.* 2019;8:883-893.
     <https://pubmed.ncbi.nlm.nih.gov/31652028/>
111. EFPIA MID3 Workgroup. **Good Practices in Model-Informed Drug Discovery and
     Development: Practice, Application, and Documentation.** *CPT Pharmacometrics Syst
     Pharmacol.* 2016;5:93-122. <https://pubmed.ncbi.nlm.nih.gov/27069774/>
112. Ribba B, Grimm HP, Agoram B, et al. **Methodologies for Quantitative Systems
     Pharmacology (QSP) Models: Design and Estimation.** *CPT Pharmacometrics Syst
     Pharmacol.* 2017;6:496-498. <https://pubmed.ncbi.nlm.nih.gov/28585415/>

**▶ Model link** — ref 107 is the basis of `CLNIC` · `VCNIC` · `QNIC` · `VPNIC` · `KE0NIC`, and
`EC5NIC` was calibrated so that the 5-15 mg/h clinical dose range lies on the **linear part** of the Emax curve
(in the draft, with `EC5NIC = 0.035 mg/L`, saturation was already reached at 3 mg/h and the guideline · intensive ·
overshoot arms became pharmacologically identical). Ref 108 is the ultra-short-acting ester
hydrolysis of clevidipine, and ref 109 the two-compartment PK of TXA and its antifibrinolytic EC50.
Refs 110-112 are mrgsolve implementation conventions and principles for QSP model documentation · validation.

---

## Appendix: what this model deliberately does not reproduce

The credibility of a QSP model shows not only in what it reproduces but also in **what it has decided not to reproduce**.

| Item | Reason |
|---|---|
| A large mRS improvement from blood pressure lowering | A result INTERACT2 · ATACH-2 did not obtain. The model also comes out small. |
| Benefit from rFVIIa · platelet transfusion | FAST · PATCH were neutral and harmful respectively. So is the model. |
| Recurrence · long-term secondary prevention | Outside the 90-day horizon, so present in the map only and not in the ODEs. |
| Individual patient prediction | Not fitted to patient-level data. Validation extends only to the direction · magnitude at trial level. |
| Molecular distinction between CAA and hypertensive pathology | Represented by the `LOC` · `CAA` nodes, but Aβ kinetics are not implemented. |

---

## Verification method

```bash
# Every PMID in this document had its title checked in exactly the way shown below.
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi\
?db=pubmed&id=<PMID>&retmode=json" | python3 -m json.tool
```

PMIDs written from memory have a high error rate. In the course of the actual look-up a substantial number of the initial
candidate PMIDs were **found to point to a different paper and were replaced**, and references that could not be
confirmed by look-up were excluded from this list.
