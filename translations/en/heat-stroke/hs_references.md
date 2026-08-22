# Heat Stroke QSP model — references

**Heat stroke — exertional (EHS) and classic / non-exertional (NEHS)**

This document collects the sources of every structural assumption and parameter that went into `hs_qsp_model.dot`, `hs_mrgsolve_model.R` and `hs_core.py`.
At the end of each section, **what this model took from the
literature concerned** is stated explicitly.

> A citation points to where the evidence lies; it does not warrant that the model is valid.
> This model is for educational and research purposes and must not be used for clinical decision-making.

---

## 1. Definitions, epidemiology and the two heat strokes

1. Bouchama A, Knochel JP. **Heat stroke.** N Engl J Med. 2002;346(25):1978-88.
   <https://pubmed.ncbi.nlm.nih.gov/12075060/>
2. Epstein Y, Yanovich R. **Heatstroke.** N Engl J Med. 2019;380(25):2449-2459.
   <https://pubmed.ncbi.nlm.nih.gov/31216400/>
3. Bouchama A, Abuyassin B, Lehe C, et al. **Classic and exertional heatstroke.**
   Nat Rev Dis Primers. 2022;8(1):8. <https://pubmed.ncbi.nlm.nih.gov/35115565/>
4. Leon LR, Bouchama A. **Heat stroke.** Compr Physiol. 2015;5(2):611-47.
   <https://pubmed.ncbi.nlm.nih.gov/25880507/>
5. Sorensen C, Hess J. **Treatment and Prevention of Heat-Related Illness.**
   N Engl J Med. 2022;387(15):1404-1413.
   <https://pubmed.ncbi.nlm.nih.gov/36239647/>
6. Argaud L, Ferry T, Le QH, et al. **Short- and long-term outcomes of heatstroke
   following the 2003 heat wave in Lyon, France.** Arch Intern Med.
   2007;167(20):2177-83. <https://pubmed.ncbi.nlm.nih.gov/17698677/>
7. Wallace RF, Kriebel D, Punnett L, et al. **Prior heat illness hospitalization
   and risk of early death.** Environ Res. 2007;104(2):290-5.
   <https://pubmed.ncbi.nlm.nih.gov/17335799/>
8. DeMartini JK, Casa DJ, Stearns R, et al. **Effectiveness of cold water
   immersion in the treatment of exertional heat stroke at the Falmouth Road
   Race.** Med Sci Sports Exerc. 2015;47(2):240-5.
   <https://pubmed.ncbi.nlm.nih.gov/24983342/>

> **In the model:** the two phenotypes (EHS/NEHS) are implemented not as separate
> diseases but as **two cases of the same heat-balance equation with
> different numerators**. The difference in mortality (EHS cooled immediately, near 0% — the Falmouth cohort, 274/274 surviving — versus NEHS
> 30–60%) is derived not from intrinsic lethality but from the clock speed of §3 and the dose of §4.

---

## 2. Thermoregulatory physiology and the heat-balance equation

9. Gagge AP, Stolwijk JA, Nishi Y. **An effective temperature scale based on a
   simple model of human physiological regulatory response.** ASHRAE Trans.
   1971;77:247-262.
10. Havenith G, Fiala D. **Thermal indices and thermophysiological modeling for
    heat stress.** Compr Physiol. 2015;6(1):255-302.
    <https://pubmed.ncbi.nlm.nih.gov/26756633/>
11. Cramer MN, Jay O. **Biophysical aspects of human thermoregulation during
    heat stress.** Auton Neurosci. 2016;196:3-13.
    <https://pubmed.ncbi.nlm.nih.gov/26971392/>
12. Cramer MN, Gagnon D, Laitano O, Crandall CG. **Human temperature regulation
    under heat stress in health, disease, and injury.** Physiol Rev.
    2022;102(4):1907-1989. <https://pubmed.ncbi.nlm.nih.gov/35679471/>
13. Kenny GP, Jay O. **Thermometry, calorimetry, and mean body temperature
    during heat stress.** Compr Physiol. 2013;3(4):1689-719.
    <https://pubmed.ncbi.nlm.nih.gov/24265242/>
14. Charkoudian N. **Skin blood flow in adult human thermoregulation.** Mayo
    Clin Proc. 2003;78(5):603-12. <https://pubmed.ncbi.nlm.nih.gov/12744548/>
15. Rowell LB. **Human cardiovascular adjustments to exercise and thermal
    stress.** Physiol Rev. 1974;54(1):75-159.
    <https://pubmed.ncbi.nlm.nih.gov/4587247/>
16. Buck AL. **New equations for computing vapor pressure and enhancement
    factor.** J Appl Meteorol. 1981;20:1527-1532.
17. Stull R. **Wet-bulb temperature from relative humidity and air
    temperature.** J Appl Meteorol Climatol. 2011;50:2267-2269.

> **In the model:** a three-node (core · muscle · skin) Gagge-family heat balance. Dry heat exchange
> `Q_dry = A(T_sk−T_o)/(R_cl+1/f_cl·h)`, evaporative capacity
> `E_max = A(P_sk−P_a)/(R_e,cl+1/f_cl·h_e)`, Lewis relation `h_e = 16.5·h_c`.
> Skin blood flow 0.3→7.5 L/min (ref 14), splanchnic blood-flow redistribution (ref 15).
> Saturation vapour pressure by Buck(16), wet-bulb temperature by the Stull(17) approximation.

---

## 3. The compensable/uncompensable boundary and the critical wet-bulb temperature (the fixed point)

18. Lind AR. **A physiological criterion for setting thermal environmental
    limits for everyday work.** J Appl Physiol. 1963;18:51-6.
    <https://pubmed.ncbi.nlm.nih.gov/13930723/>
19. Belding HS, Hatch TF. **Index for evaluating heat stress in terms of
    resulting physiological strains.** Heat Pip Air Cond. 1955;27:129-136.
20. Sherwood SC, Huber M. **An adaptability limit to climate change due to heat
    stress.** Proc Natl Acad Sci USA. 2010;107(21):9552-5.
    <https://pubmed.ncbi.nlm.nih.gov/20439769/>
21. Vecellio DJ, Wolf ST, Cottle RM, Kenney WL. **Evaluating the 35°C wet-bulb
    temperature adaptability threshold for young, healthy subjects
    (PSU HEAT Project).** J Appl Physiol. 2022;132(2):340-345.
    <https://pubmed.ncbi.nlm.nih.gov/34913738/>
22. Wolf ST, Kenney WL. **The dry heat is hotter: the critical environmental
    limits of young adults in humid vs dry heat.** Am J Physiol Regul Integr
    Comp Physiol. 2023;324(5):R601-R608.
    <https://pubmed.ncbi.nlm.nih.gov/36939211/>
23. Foster J, Smallcombe JW, Hodder S, et al. **An advanced empirical model for
    quantifying the impact of heat and climate change on human physical work
    capacity.** Int J Biometeorol. 2021;65(7):1215-1229.
    <https://pubmed.ncbi.nlm.nih.gov/33674931/>
24. Vanos J, Guzman-Echavarria G, Baldwin JW, et al. **A physiological approach
    for assessing human survivability and liveability to heat in a changing
    climate.** Nat Commun. 2023;14(1):7653.
    <https://pubmed.ncbi.nlm.nih.gov/37993427/>

> **In the model:** this model's **result no. 1**. Compensability is defined not as a problem of simulation but
> as the question of **whether a fixed point exists in the heat-balance equation**.
> The calculation gives a critical wet-bulb temperature at rest of **35.7 °C** — the 35 °C
> limit of Sherwood-Huber(20) falls out of the equation even though it was never entered as a parameter — and at 900 W of work it
> drops to **about 21 °C**. That reproduces the measurements of Vecellio(21) and Wolf(22): "the
> actual critical wet-bulb temperature of young healthy adults is far below 35 °C".
> **A refutation obtained along the way:** at every condition above this boundary `E_actual = E_max,env`
> holds, and sweating capacity (unacclimatised 668 W / acclimatised 1233 W) is not once the limiting
> factor. So in this model **acclimatisation cannot move the boundary** (mean shift
> −0.35 °C). It becomes a testable structural claim that the gain from acclimatisation lies not in the boundary but in the transient and in
> cardiovascular reserve.

---

## 4. Thermal dose: CEM43 (Sapareto-Dewey)

25. Sapareto SA, Dewey WC. **Thermal dose determination in cancer therapy.**
    Int J Radiat Oncol Biol Phys. 1984;10(6):787-800.
    <https://pubmed.ncbi.nlm.nih.gov/6547421/>
26. Dewhirst MW, Viglianti BL, Lora-Michiels M, et al. **Basic principles of
    thermal dosimetry and thermal thresholds for tissue damage from
    hyperthermia.** Int J Hyperthermia. 2003;19(3):267-94.
    <https://pubmed.ncbi.nlm.nih.gov/12745972/>
27. van Rhoon GC, Samaras T, Yarmolenko PS, et al. **CEM43°C thermal dose
    thresholds: a potential guide for magnetic resonance radiofrequency
    exposure levels?** Eur Radiol. 2013;23(8):2215-27.
    <https://pubmed.ncbi.nlm.nih.gov/23553588/>
28. Yarmolenko PS, Moon EJ, Landon C, et al. **Thresholds for thermal damage to
    normal tissues: an update.** Int J Hyperthermia. 2011;27(4):320-43.
    <https://pubmed.ncbi.nlm.nih.gov/21591897/>
29. Wright NT. **On a relationship between the Arrhenius parameters from
    thermal damage studies.** J Biomech Eng. 2003;125(2):300-4.
    <https://pubmed.ncbi.nlm.nih.gov/12751294/>

> **In the model:** this model's **result no. 2**. The injury variable is taken not as peak body temperature but as
> the **cumulative thermal dose CEM43 = ∫R^(43−T_c)dt** (R = 0.25 below 43 °C, 0.50 above).
> That single line forces "time and temperature are not exchanged 1:1" —
> 0.016/min at 40 °C, 0.250/min at 42 °C, 1.000/min at 43 °C.
> Applying to heat stroke this dose law, established in tumour hyperthermia, is the core
> transplant of this model, and it is what makes it possible to place cooling delay and
> cooling modality **on a single axis**.

---

## 5. Cooling: modality and rate

30. Casa DJ, McDermott BP, Lee EC, et al. **Cold water immersion: the gold
    standard for exertional heatstroke treatment.** Exerc Sport Sci Rev.
    2007;35(3):141-9. <https://pubmed.ncbi.nlm.nih.gov/17620933/>
31. Casa DJ, Armstrong LE, Kenny GP, et al. **Exertional heat stroke: new
    concepts regarding cause and care.** Curr Sports Med Rep. 2012;11(3):115-23.
    <https://pubmed.ncbi.nlm.nih.gov/22580488/>
32. McDermott BP, Casa DJ, Ganio MS, et al. **Acute whole-body cooling for
    exercise-induced hyperthermia: a systematic review.** J Athl Train.
    2009;44(1):84-93. <https://pubmed.ncbi.nlm.nih.gov/19180223/>
33. Zhang Y, Davis JK, Casa DJ, Bishop PA. **Optimizing cold water immersion for
    exercise-induced hyperthermia: a meta-analysis.** Med Sci Sports Exerc.
    2015;47(11):2464-72. <https://pubmed.ncbi.nlm.nih.gov/25848900/>
34. Proulx CI, Ducharme MB, Kenny GP. **Effect of water temperature on cooling
    efficiency during hyperthermia in humans.** J Appl Physiol.
    2003;94(4):1317-23. <https://pubmed.ncbi.nlm.nih.gov/12626467/>
35. Luhring KE, Butts CL, Smith CR, et al. **Cooling effectiveness of a modified
    cold-water immersion method after exercise-induced hyperthermia.**
    J Athl Train. 2016;51(11):946-951.
    <https://pubmed.ncbi.nlm.nih.gov/27834504/>
36. Hosokawa Y, Adams WM, Belval LN, et al. **Tarp-assisted cooling as a method
    of whole-body cooling in hyperthermic individuals.** Ann Emerg Med.
    2017;69(3):347-352. <https://pubmed.ncbi.nlm.nih.gov/27522596/>
37. Butts CL, McDermott BP, Buening BJ, et al. **Physiologic and perceptual
    responses to cold-shower cooling after exercise-induced hyperthermia.**
    J Athl Train. 2016;51(3):252-7.
    <https://pubmed.ncbi.nlm.nih.gov/26967831/>
38. Gaudio FG, Grissom CK. **Cooling methods in heat stroke.** J Emerg Med.
    2016;50(4):607-16. <https://pubmed.ncbi.nlm.nih.gov/26525947/>
39. Filep EM, Murata Y, Endres BD, et al. **Exertional heat stroke, modality
    cooling rate, and survival outcomes: a systematic review.** Medicina
    (Kaunas). 2020;56(11):589. <https://pubmed.ncbi.nlm.nih.gov/33167418/>
40. Weiner JS, Khogali M. **A physiological body-cooling unit for treatment of
    heat stroke.** Lancet. 1980;1(8167):507-9.
    <https://pubmed.ncbi.nlm.nih.gov/6102236/>
41. Hoedemaekers CW, Ezzahti M, Gerritsen A, van der Hoeven JG. **Comparison of
    cooling methods to induce and maintain normo- and hypothermia in ICU
    patients.** Crit Care. 2007;11(4):R91.
    <https://pubmed.ncbi.nlm.nih.gov/17718920/>
42. Kim F, Olsufka M, Longstreth WT Jr, et al. **Pilot randomized clinical trial
    of prehospital induction of mild hypothermia with cold normal saline.**
    Circulation. 2007;115(24):3064-70.
    <https://pubmed.ncbi.nlm.nih.gov/17548731/>
43. Bernard S, Buist M, Monteiro O, Smith K. **Induced hypothermia using large
    volume, ice-cold intravenous fluid in comatose survivors of out-of-hospital
    cardiac arrest.** Resuscitation. 2003;56(1):9-13.
    <https://pubmed.ncbi.nlm.nih.gov/12505732/>
44. Belval LN, Casa DJ, Adams WM, et al. **Consensus statement — prehospital
    care of exertional heat stroke.** Prehosp Emerg Care. 2018;22(3):392-397.
    <https://pubmed.ncbi.nlm.nih.gov/29336710/>

> **In the model:** each cooling modality is expressed as **a single lumped conductance UA (W/K)**, and was
> fitted numerically in `hs_calibrate.py` against the published core cooling rates of the
> literature above: ice-water immersion 2 °C → UA 23.5 (0.22 °C/min), cold-water immersion 14 °C → 27.7
> (0.17), tarp-assisted cooling → 20.1 (0.14), cold shower → 15.9 (0.10), evaporative-convective →
> 16.2 (0.08), intravascular catheter → 4.6 (0.06), ice packs → 1.6 (0.032),
> passive → 0.016. The point is that the modality has been turned into **a single potency number** —
> it puts "which cooler" and "when it is started" on the same axis.
> Two litres of 4 °C crystalloid is implemented not as a conductance but as **a one-off enthalpy sink**
> (perfusion-weighted distribution), and the calculated core drop of 1.4–1.7 °C agrees with the
> observations of Kim(42) and Bernard(43).
> **Where the model over-predicts:** it predicts the ratio of the cooling rates of 2 °C versus 14 °C water as 1.55,
> whereas the meta-analytic measurement is 1.16–1.29. The cause is that the conductance is linear and that shivering
> heat production offset and incomplete immersion were not modelled;
> it is flagged in `hs_verification_output.txt` as this model's most exposed thermodynamic
> prediction. Because the scenarios use a separately fitted UA for each modality, the reported
> cooling rates themselves do agree with the observations.

---

## 6. Cooling delay determines the outcome (time to cooling)

45. Heled Y, Rav-Acha M, Shani Y, et al. **The "golden hour" for heatstroke
    treatment.** Mil Med. 2004;169(3):184-6.
    <https://pubmed.ncbi.nlm.nih.gov/15080235/>
46. Vicario SJ, Okabajue R, Haltom T. **Rapid cooling in classic heatstroke:
    effect on mortality rates.** Am J Emerg Med. 1986;4(5):394-8.
    <https://pubmed.ncbi.nlm.nih.gov/3741548/>
47. Demartini JK, Casa DJ, Belval LN, et al. **Environmental conditions and the
    occurrence of exertional heat illnesses and exertional heat stroke at the
    Falmouth Road Race.** J Athl Train. 2014;49(4):478-85.
    <https://pubmed.ncbi.nlm.nih.gov/24960519/>
48. Rublee C, Dresser C, Giudice C, et al. **Evidence-based heatstroke
    management in the emergency department.** West J Emerg Med.
    2021;22(2):186-195. <https://pubmed.ncbi.nlm.nih.gov/33856299/>
49. Pryor RR, Roth RN, Suyama J, Hostler D. **Exertional heat illness:
    emerging concepts and advances in prehospital care.** Prehosp Disaster Med.
    2015;30(3):297-305. <https://pubmed.ncbi.nlm.nih.gov/25868636/>

> **In the model:** this model's **result no. 3 — the exchange rate**.
> In a patient who collapsed at 42.0 °C, **upgrading the modality** from evaporative cooling to ice-water immersion is
> worth **4.6–6.0 minutes of delay** on the CEM43 scale, and this exchange rate is almost independent of the
> collapse temperature (41/42/43 °C). The absolute risk, by contrast, grows exponentially.
> This is exactly the quantitative content of "cool first, transport second"(44, 48):
> modality is worth a factor of a few, delay has no upper bound.

---

## 7. Gastrointestinal barrier, endotoxaemia and splanchnic ischaemia

50. Lambert GP. **Stress-induced gastrointestinal barrier dysfunction and its
    inflammatory effects.** J Anim Sci. 2009;87(14 Suppl):E101-8.
    <https://pubmed.ncbi.nlm.nih.gov/18791134/>
51. Lim CL, Mackinnon LT. **The roles of exercise-induced immune system
    disturbances in the pathology of heat stroke.** Sports Med.
    2006;36(1):39-64. <https://pubmed.ncbi.nlm.nih.gov/16445310/>
52. Hall DM, Buettner GR, Oberley LW, et al. **Mechanisms of circulatory and
    intestinal barrier dysfunction during whole body hyperthermia.** Am J
    Physiol Heart Circ Physiol. 2001;280(2):H509-21.
    <https://pubmed.ncbi.nlm.nih.gov/11158946/>
53. Dokladny K, Zuhl MN, Moseley PL. **Intestinal epithelial barrier function
    and tight junction proteins with heat and exercise.** J Appl Physiol.
    2016;120(6):692-701. <https://pubmed.ncbi.nlm.nih.gov/26359485/>
54. Ogden HB, Child RB, Fallowfield JL, et al. **The gastrointestinal exertional
    heat stroke paradigm: pathophysiology, assessment, severity, aetiology and
    nutritional countermeasures.** Nutrients. 2020;12(2):537.
    <https://pubmed.ncbi.nlm.nih.gov/32093001/>
55. Bouchama A, Parhar RS, el-Yazigi A, et al. **Endotoxemia and release of
    tumor necrosis factor and interleukin 1 alpha in acute heatstroke.**
    J Appl Physiol. 1991;70(6):2640-4.
    <https://pubmed.ncbi.nlm.nih.gov/1885459/>
56. Snipe RMJ, Khoo A, Kitic CM, et al. **The impact of exertional-heat stress
    on gastrointestinal integrity, systemic endotoxin and cytokine profile.**
    Eur J Appl Physiol. 2018;118(2):389-400.
    <https://pubmed.ncbi.nlm.nih.gov/29234915/>

> **In the model:** the **cutaneous steal** by which skin vasodilatation (the very mechanism that dumps the heat)
> takes splanchnic blood flow from 1.5 → 0.3 L/min is implemented explicitly.
> Enterocyte injury is the sum of two terms, ischaemia (ISCH²) and direct thermal injury (CEM43), and the event that opens
> one door **closes the exit at the same time** — because the endotoxin clearance of Kupffer cells is
> proportional both to liver function and to splanchnic perfusion.

---

## 8. Cytokines, HMGB1 and the commitment switch

57. Bouchama A, Hammami MM, Haq A, et al. **Evidence for endothelial cell
    activation/injury in heatstroke.** Crit Care Med. 1996;24(7):1173-8.
    <https://pubmed.ncbi.nlm.nih.gov/8674331/>
58. Bouchama A, Roberts G, Al Mohanna F, et al. **Inflammatory, hemostatic, and
    clinical changes in a baboon experimental model for heatstroke.**
    J Appl Physiol. 2005;98(2):697-705.
    <https://pubmed.ncbi.nlm.nih.gov/15475605/>
59. Leon LR, Helwig BG. **Heat stroke: role of the systemic inflammatory
    response.** J Appl Physiol. 2010;109(6):1980-8.
    <https://pubmed.ncbi.nlm.nih.gov/20522730/>
60. Tong H, Tang Y, Chen Y, et al. **HMGB1 activity inhibition alleviating liver
    injury in heatstroke.** J Trauma Acute Care Surg. 2013;74(3):801-7.
    <https://pubmed.ncbi.nlm.nih.gov/23425739/>
61. Tong HS, Tang YQ, Chen Y, et al. **Early elevated HMGB1 level predicting the
    outcome in exertional heatstroke.** J Trauma. 2011;71(4):808-14.
    <https://pubmed.ncbi.nlm.nih.gov/21460740/>
62. Andersson U, Tracey KJ. **HMGB1 is a therapeutic target for sterile
    inflammation and infection.** Annu Rev Immunol. 2011;29:139-62.
    <https://pubmed.ncbi.nlm.nih.gov/21219181/>
63. Scaffidi P, Misteli T, Bianchi ME. **Release of chromatin protein HMGB1 by
    necrotic cells triggers inflammation.** Nature. 2002;418(6894):191-5.
    <https://pubmed.ncbi.nlm.nih.gov/12110890/>
64. Dehbi M, Baturcam E, Eldali A, et al. **Hsp-72, a candidate prognostic
    indicator of heatstroke.** Cell Stress Chaperones. 2010;15(5):593-603.
    <https://pubmed.ncbi.nlm.nih.gov/20205026/>
65. Ferat-Osorio E, Sánchez-Anaya A, Gutiérrez-Mendoza M, et al. **Heat shock
    protein 70 down-regulates the production of TLR-mediated inflammatory
    cytokines.** J Inflamm (Lond). 2014;11:19.
    <https://pubmed.ncbi.nlm.nih.gov/25053921/>

> **In the model:** this model's **result no. 4 — regime 3**. HMGB1 is taken as the commitment variable and
> implemented as a **saddle-node bistability**:
> `dH/dt = drive + A·H³/(H³+K³) − c_eff·H`.
> The fixed points are OFF = 0, **unstable = 10.91 ng/mL**, ON = 50.3 ng/mL.
> Necrotic cells release HMGB1(63) and HMGB1 in turn creates further necrosis
> through RAGE/TLR4(60, 62); when that positive feedback is joined by the two conditions that production saturates
> and clearance is linear, two stable states necessarily arise.
> **The key calculation:** the **dose that produces commitment is 7.8 ± 0.2 CEM43, almost independent of
> the cooling modality** (spread 0.69), whereas the **delay that produces commitment ranges from 17 minutes (ice packs)
> to 42 minutes (ice-water immersion)**. The dose is a property of the patient; the delay that can be
> tolerated is a property of the cooler.
> HSP70 enters as a protective term dividing the thermal dose(64, 65).

---

## 9. Coagulopathy and DIC

66. Bouchama A, Bridey F, Hammami MM, et al. **Activation of coagulation and
    fibrinolysis in heatstroke.** Thromb Haemost. 1996;76(6):909-15.
    <https://pubmed.ncbi.nlm.nih.gov/8972008/>
67. Hifumi T, Kondo Y, Shimizu K, Miyake Y. **Heat stroke.** J Intensive Care.
    2018;6:30. <https://pubmed.ncbi.nlm.nih.gov/29850022/>
68. Iba T, Levy JH, Warkentin TE, et al. **Diagnosis and management of
    sepsis-induced coagulopathy and disseminated intravascular coagulation.**
    J Thromb Haemost. 2019;17(11):1989-1994.
    <https://pubmed.ncbi.nlm.nih.gov/31410983/>
69. Taylor FB Jr, Toh CH, Hoots WK, et al. **Towards definition, clinical and
    laboratory criteria, and a scoring system for disseminated intravascular
    coagulation (ISTH).** Thromb Haemost. 2001;86(5):1327-30.
    <https://pubmed.ncbi.nlm.nih.gov/11816725/>
70. Ogura Y, Sakurai R, Kawakami R, et al. **Heat-stroke-induced DIC and
    recombinant thrombomodulin.** Acute Med Surg. 2018;5(4):327-334.
    <https://pubmed.ncbi.nlm.nih.gov/30338078/>
71. Yokobori S, Koido Y, Shishido H, et al. **Feasibility and safety of
    intravascular temperature management for severe heat stroke (ICE-HEAT).**
    Crit Care Med. 2018;46(7):e670-e676.
    <https://pubmed.ncbi.nlm.nih.gov/29652723/>
72. Johansson PI, Stensballe J, Ostrowski SR. **Shock-induced endotheliopathy
    (SHINE) in acute critical illness.** Crit Care. 2017;21(1):25.
    <https://pubmed.ncbi.nlm.nih.gov/28173843/>

> **In the model:** implemented as a consumptive coagulopathy — thrombin generation → consumption of
> fibrinogen · platelets · protein C → rising D-dimer. **The ISTH overt-DIC score and SOFA · GCS are all
> computed as outputs only, not as state variables** — integrating a score would let the score
> drive the pathophysiology, which is a circular argument.
> Fibrinogen is the only factor whose consumption term and acute-phase synthesis term have **opposite signs**.

---

## 10. Drugs: antipyretics, dantrolene, steroids, thrombomodulin

73. Bouchama A, Cafege A, Devol EB, et al. **Ineffectiveness of dantrolene
    sodium in the treatment of heatstroke.** Crit Care Med. 1991;19(2):176-80.
    <https://pubmed.ncbi.nlm.nih.gov/1989755/>
74. Krause T, Gerbershagen MU, Fiege M, et al. **Dantrolene — a review of its
    pharmacology, therapeutic use and new developments.** Anaesthesia.
    2004;59(4):364-73. <https://pubmed.ncbi.nlm.nih.gov/15023108/>
75. Rosenberg H, Pollock N, Schiemann A, et al. **Malignant hyperthermia: a
    review.** Orphanet J Rare Dis. 2015;10:93.
    <https://pubmed.ncbi.nlm.nih.gov/26238698/>
76. Prescott LF. **Paracetamol, alcohol and the liver.** Br J Clin Pharmacol.
    2000;49(4):291-301. <https://pubmed.ncbi.nlm.nih.gov/10759684/>
77. Mitchell JR, Jollow DJ, Potter WZ, et al. **Acetaminophen-induced hepatic
    necrosis. IV. Protective role of glutathione.** J Pharmacol Exp Ther.
    1973;187(1):211-7. <https://pubmed.ncbi.nlm.nih.gov/4746329/>
78. Whitcomb DC, Block GD. **Association of acetaminophen hepatotoxicity with
    fasting and ethanol use.** JAMA. 1994;272(23):1845-50.
    <https://pubmed.ncbi.nlm.nih.gov/7990219/>
79. Saito T, Maruyama I, Shimazaki S, et al. **Efficacy and safety of
    recombinant human soluble thrombomodulin (ART-123) in disseminated
    intravascular coagulation.** J Thromb Haemost. 2007;5(1):31-41.
    <https://pubmed.ncbi.nlm.nih.gov/17059423/>
80. Abeyama K, Stern DM, Ito Y, et al. **The N-terminal domain of
    thrombomodulin sequesters high-mobility group-B1 protein: a novel
    antiinflammatory mechanism.** J Clin Invest. 2005;115(5):1267-74.
    <https://pubmed.ncbi.nlm.nih.gov/15841214/>
81. Vaity C, Al-Subaie N, Cecconi M. **Cooling techniques for hyperthermia in
    the intensive care unit.** Crit Care. 2015;19:103.
    <https://pubmed.ncbi.nlm.nih.gov/25886948/>
82. Niven DJ, Stelfox HT, Laupland KB. **Antipyretic therapy in febrile
    critically ill adults: a systematic review and meta-analysis.**
    J Crit Care. 2013;28(3):303-10.
    <https://pubmed.ncbi.nlm.nih.gov/23159136/>
83. Sund-Levander M, Grodzinsky E. **Assessment of body temperature measurement
    options.** Br J Nurs. 2013;22(15):880-8.
    <https://pubmed.ncbi.nlm.nih.gov/24006812/>

> **In the model:** this model's **result no. 5 — the antipyretic argument**. Fever raises the **set point**
> and the body defends it, but in heat stroke the set point stays at 37 °C and the effectors are
> already saturated. Antipyretics are therefore implemented as acting **on the threshold only**, and
> the set-point sensitivity of dTc/dt is large in the compensable region and **collapses to zero in
> the uncompensable region** — an antipyretic works only in the region that is not heat stroke(82).
> The cost is billed separately: NAPQI in a liver where CYP2E1 is induced and glutathione is
> already depleted(76-78), and the renal prostaglandin blockade and platelet-function
> inhibition of NSAIDs.
> Dantrolene acts through the RyR1 mechanism(74, 75) — that is, **the mechanism of malignant hyperthermia, not the
> mechanism of this one** — so it cuts only ≤12% of muscle heat production, reproducing the null
> observed in Bouchama's RCT(73).
> Thrombomodulin alfa is, in this model, **the only drug that touches the commitment variable**:
> because its lectin-like domain degrades HMGB1(80) it adds to HMGB1 clearance, and on the
> saddle-node calculation **beyond +0.00192/min the ON state itself disappears**.
> Its therapeutic window is therefore set by the switch, not by the thermometer — a testable
> structural prediction(70, 79).

---

## 11. Acclimatisation, dehydration and risk factors

84. Périard JD, Racinais S, Sawka MN. **Adaptations and mechanisms of human heat
    acclimation.** Scand J Med Sci Sports. 2015;25 Suppl 1:20-38.
    <https://pubmed.ncbi.nlm.nih.gov/25943654/>
85. Périard JD, Eijsvogels TMH, Daanen HAM. **Exercise under heat stress:
    thermoregulation, hydration, performance implications, and mitigation
    strategies.** Physiol Rev. 2021;101(4):1873-1979.
    <https://pubmed.ncbi.nlm.nih.gov/33829868/>
86. Sawka MN, Burke LM, Eichner ER, et al. **ACSM position stand: exercise and
    fluid replacement.** Med Sci Sports Exerc. 2007;39(2):377-90.
    <https://pubmed.ncbi.nlm.nih.gov/17277604/>
87. Cheuvront SN, Kenefick RW. **Dehydration: physiology, assessment, and
    performance effects.** Compr Physiol. 2014;4(1):257-85.
    <https://pubmed.ncbi.nlm.nih.gov/24692140/>
88. Kenney WL, Craighead DH, Alexander LM. **Heat waves, aging, and human
    cardiovascular health.** Med Sci Sports Exerc. 2014;46(10):1891-9.
    <https://pubmed.ncbi.nlm.nih.gov/24598696/>
89. Semenza JC, Rubin CH, Falter KH, et al. **Heat-related deaths during the
    July 1995 heat wave in Chicago.** N Engl J Med. 1996;335(2):84-90.
    <https://pubmed.ncbi.nlm.nih.gov/8649494/>
90. Westaway K, Frank O, Husband A, et al. **Medicines can affect
    thermoregulation and accentuate the risk of dehydration and heat-related
    illness during hot weather.** J Clin Pharm Ther. 2015;40(4):363-7.
    <https://pubmed.ncbi.nlm.nih.gov/25917210/>
91. Bongers CCWG, Alsady M, Nijenhuis T, et al. **Impact of acute versus
    prolonged exercise and dehydration on kidney function and injury.**
    Physiol Rep. 2018;6(11):e13734. <https://pubmed.ncbi.nlm.nih.gov/29890037/>

> **In the model:** acclimatisation is tied to a single knob (`ACCLIM`) as SWMAX ×1.85, sweating threshold −0.35 °C,
> plasma volume +11.5% and a raised HSP baseline(84, 85).
> Dehydration cuts sweating by 5.5% per 1% of body-weight deficit(86, 87) and carries through to
> plasma volume · MAP · GFR. Anticholinergic burden(90) enters as `FSW_DRUG`, ageing as `FSW_AGE` · `FVD_AGE`
> — exactly the risk-factor structure identified in the Chicago 1995 heat wave(89).
> In the model **a single anticholinergic flips whether a fixed point exists**: in the same 40 °C/55% RH
> room, an older person without an anticholinergic creeps slowly upwards above the boundary (38.1 °C at
> 72 hours), while one taking it loses the fixed point and runs away.

---

## 12. Organ injury, rhabdomyolysis and prognosis

92. Argaud L, Ferry T, Le QH, et al. (see ref 6 — long-term outcome)
93. Pease S, Bouadma L, Kermarrec N, et al. **Early organ dysfunction course,
    cooling time and outcome in classic heatstroke.** Intensive Care Med.
    2009;35(8):1454-8. <https://pubmed.ncbi.nlm.nih.gov/19404609/>
94. Varghese GM, John G, Thomas K, et al. **Predictors of multi-organ
    dysfunction in heatstroke.** Emerg Med J. 2005;22(3):185-7.
    <https://pubmed.ncbi.nlm.nih.gov/15735266/>
95. Hifumi T, Kondo Y, Shimazaki J, et al. **Prognostic significance of
    disseminated intravascular coagulation in patients with heat stroke.**
    J Intensive Care. 2018;6:33. <https://pubmed.ncbi.nlm.nih.gov/29942522/>
96. Garcin JM, Bronstein JA, Cremades S, et al. **Acute liver failure is
    frequent during heat stroke.** World J Gastroenterol. 2008;14(1):158-9.
    <https://pubmed.ncbi.nlm.nih.gov/18176982/>
97. Zeller L, Novack V, Barski L, et al. **Exertional heatstroke: clinical
    characteristics, diagnostic and therapeutic considerations.** Eur J Intern
    Med. 2011;22(3):296-9. <https://pubmed.ncbi.nlm.nih.gov/21570652/>
98. Bosch X, Poch E, Grau JM. **Rhabdomyolysis and acute kidney injury.**
    N Engl J Med. 2009;361(1):62-72.
    <https://pubmed.ncbi.nlm.nih.gov/19571284/>
99. Lawton EM, Pearce H, Gabb GM. **Review article: environmental heatstroke and
    long-term clinical neurological outcomes.** Emerg Med Australas.
    2019;31(2):163-173. <https://pubmed.ncbi.nlm.nih.gov/30632286/>
100. Yaqub BA. **Neurologic manifestations of heatstroke at the Mecca
     pilgrimage.** Neurology. 1987;37(6):1004-6.
     <https://pubmed.ncbi.nlm.nih.gov/3587615/>

> **In the model:** liver injury is the sum of three terms, CEM43 · ischaemia · NAPQI, and the half-lives are set so that AST/ALT
> peak at 48–72 hours(96). Rhabdomyolysis is driven by **the CEM43 of the muscle node,
> not of the core** — during exercise the muscle is 1–2 °C above what a rectal probe
> reads, so the reason CK is far higher in EHS than in NEHS at the same core temperature falls out
> automatically. The observation that cerebellar Purkinje cells are the most heat-sensitive
> neurons(99, 100) is reflected as a low threshold in the CNS injury term.

---

## 13. Clinical guidelines and consensus statements

101. Casa DJ, DeMartini JK, Bergeron MF, et al. **National Athletic Trainers'
     Association Position Statement: Exertional Heat Illnesses.** J Athl Train.
     2015;50(9):986-1000. <https://pubmed.ncbi.nlm.nih.gov/26381473/>
102. Roberts WO, Armstrong LE, Sawka MN, et al. **ACSM Expert Consensus
     Statement on Exertional Heat Illness: Recognition, Management, and Return
     to Activity.** Curr Sports Med Rep. 2021;20(9):470-484.
     <https://pubmed.ncbi.nlm.nih.gov/34524191/>
103. Lipman GS, Gaudio FG, Eifling KP, et al. **Wilderness Medical Society
     Clinical Practice Guidelines for the Prevention and Treatment of
     Heat Illness: 2019 Update.** Wilderness Environ Med. 2019;30(4S):S33-S46.
     <https://pubmed.ncbi.nlm.nih.gov/31221601/>
104. Epstein Y, Roberts WO. **The pathophysiology of heat stroke: an integrative
     view of the final common pathway.** Scand J Med Sci Sports.
     2011;21(6):742-8. <https://pubmed.ncbi.nlm.nih.gov/21631615/>
105. Hifumi T, Kondo Y, Shimazaki J, et al. **Japanese Association for Acute
     Medicine Heatstroke and Hypothermia Surveillance Committee: heatstroke
     STUDY.** Acute Med Surg. 2019;6(3):270-278.
     <https://pubmed.ncbi.nlm.nih.gov/31304029/>
106. Périard JD, DeGroot D, Jay O. **Exertional heat stroke in sport and the
     military: epidemiology and mitigation.** Exp Physiol.
     2022;107(10):1111-1121. <https://pubmed.ncbi.nlm.nih.gov/35953078/>

> **In the model:** the cooling stop threshold of 38.6 °C(101), "cool first, transport
> second"(101-103), and work-rest cycling and WBGT-based management(106) are reflected in the scenarios and in the
> parameter defaults.

---

## 14. This model's exposed predictions and limitations

The points where the model **diverges** from the literature, or cannot be verified, are set out here. These are
where it will be wrong first.

| # | Prediction / limitation | Status |
|---|-----------|------|
| 1 | Ratio of the cooling rates of 2 °C versus 14 °C water **1.55** (meta-analysis 1.16–1.29) | **Over-prediction.** Linear conductance; shivering offset not modelled; incomplete immersion not considered. The scenarios use a separately fitted UA per modality, so the reported rates themselves do agree |
| 2 | Acclimatisation **cannot move the compensability boundary at all** | Counter-intuitive. Inside the model it is inevitable, since `E_actual = E_max,env`; but if WMAX 0.85 or the convection assumptions are generous, sweating could in reality become the limiting factor in dry heat |
| 3 | Commitment dose **7.8 CEM43** | Calibrated against the epidemiology (survival when cooled within 30 minutes, not when delayed by 60), not an independent measurement. Calibration of the absolute HMGB1 values relies on a small cohort(61) |
| 4 | rTM makes the ON state **disappear** | A direct consequence of the saddle-node calculation. The only clinical evidence is a retrospective registry study(70), and the SCARLET trial in sepsis DIC was null |
| 5 | Resting core equilibrium at **36.8 °C** (normal 37.0) | Because active cold defence was not included. Biased in the direction of delaying the onset of NEHS slightly |
| 6 | No inter-individual variability | A single deterministic patient. Risk factors such as genetic susceptibility (RYR1 variants), recent fever, sleep deprivation and obesity enter only as parameters |

---

## 15. QSP methodology

107. Gadkar K, Kirouac DC, Mager DE, et al. **A six-stage workflow for robust
     application of systems pharmacology.** CPT Pharmacometrics Syst Pharmacol.
     2016;5(5):235-49. <https://pubmed.ncbi.nlm.nih.gov/27299936/>
108. Baron KT, Gastonguay MR. **Simulation from ODE-based population PK/PD and
     systems pharmacology models in R with mrgsolve.** J Pharmacokinet
     Pharmacodyn. 2015;42:S84-S85.
109. Musante CJ, Ramanujan S, Schmidt BJ, et al. **Quantitative systems
     pharmacology: a case for disease models.** Clin Pharmacol Ther.
     2017;101(1):24-27. <https://pubmed.ncbi.nlm.nih.gov/27709613/>
110. Ribba B, Grimm HP, Agoram B, et al. **Methodologies for quantitative
     systems pharmacology (QSP) models.** CPT Pharmacometrics Syst Pharmacol.
     2017;6(8):496-498. <https://pubmed.ncbi.nlm.nih.gov/28643452/>

---

**110 references in total.** All links are to PubMed or PMC.
Parameter sources are also marked inline in the `$PARAM` comments of `hs_mrgsolve_model.R`
and in the `P0` dictionary of `hs_core.py`.
