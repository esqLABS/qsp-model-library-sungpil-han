# Acute organophosphorus insecticide poisoning — References
# Acute Organophosphorus (OP) Insecticide Self-Poisoning — References

This document organises, section by section, **the basis of every parameter and structural assumption**
that went into `op_qsp_model.dot` (the mechanistic map), `op_mrgsolve_model.R` (the 51-state ODE
model), `op_reference_model.py` (the Python re-implementation used for verification) and
`op_shiny_app.R` (the dashboard). The brackets at the end of each entry show which number in the
model that reference went into.

> **Suggested reading order** — the model's core arithmetic (the oxime ceiling Ω) is in §2, the
> comparison against clinical outcome in §7 · §8, and "why the randomised trial was negative" in §8 · §13.

---

## 1. Epidemiology and disease burden

Deliberate self-poisoning with organophosphorus insecticides accounts for a substantial share of
suicide deaths worldwide, and the best-established fact in the field is that its case fatality is
determined not by pharmacology but by **which formulations are on the market**. That fact is
reflected in the model's `REG` (regulation) node and in the virtual trial results of §13.

1. Mew EJ, Padmanathan P, Konradsen F, et al. The global burden of fatal self-poisoning with pesticides 2006-15: systematic review. *J Affect Disord* 2017;219:93-104. [Global mortality burden · model preamble] <https://pubmed.ncbi.nlm.nih.gov/28535450/>
2. Gunnell D, Eddleston M, Phillips MR, Konradsen F. The global distribution of fatal pesticide self-poisoning: systematic review. *BMC Public Health* 2007;7:357. [Annual mortality estimate] <https://pubmed.ncbi.nlm.nih.gov/18154668/>
3. Eddleston M, Buckley NA, Eyer P, Dawson AH. Management of acute organophosphorus pesticide poisoning. *Lancet* 2008;371:597-607. [The overall clinical frame; the priority order of atropine · oxime · ventilation] <https://pubmed.ncbi.nlm.nih.gov/17706760/>
4. Eddleston M, Eyer P, Worek F, et al. Differences between organophosphorus insecticides in human self-poisoning: a prospective cohort study. *Lancet* 2005;366:1452-9. [The dimethyl vs diethyl difference in case fatality · the basis for the `subclass` parameter] <https://pubmed.ncbi.nlm.nih.gov/16243090/>
5. Eddleston M, Eyer P, Worek F, et al. Predicting outcome using butyrylcholinesterase activity in organophosphorus pesticide self-poisoning. *QJM* 2008;101:467-74. [The prognostic value of BChE · the `B_free` output] <https://pubmed.ncbi.nlm.nih.gov/18375475/>
6. Dawson AH, Eddleston M, Senarathna L, et al. Acute human lethal toxicity of agricultural pesticides: a prospective cohort study. *PLoS Med* 2010;7:e1000357. [Case fatality ranking by product] <https://pubmed.ncbi.nlm.nih.gov/20968586/>
7. Manuweera G, Eddleston M, Egodage S, Buckley NA. Do targeted bans of insecticides to prevent deaths from self-poisoning result in reduced agricultural output? *Environ Health Perspect* 2008;116:492-5. [The feasibility of a regulatory intervention] <https://pubmed.ncbi.nlm.nih.gov/18414632/>
8. Knipe DW, Chang SS, Dawson A, et al. Suicide prevention through means restriction: impact of the 2008-2011 pesticide restrictions on suicide in Sri Lanka. *PLoS One* 2017;12:e0172893. [The country-level effect of a ban] <https://pubmed.ncbi.nlm.nih.gov/28264041/>

---

## 2. Esterase kinetics — the core of the model

The **k_i, k_a, k_s and k_r2 values entering the model's three-state switch (E → EP → EP-aged) all
come from the human AChE measurements of the Worek group** below. In particular, the fact that oxime
reactivation **saturates** in the form `k_r2 = k_r_max/(K_D + X)` is what generates the model's
ceiling `Ω/(1+Ω)`.

9. Worek F, Thiermann H, Szinicz L, Eyer P. Kinetic analysis of interactions between human acetylcholinesterase, structurally different organophosphorus compounds and oximes. *Biochem Pharmacol* 2004;68:2237-48. [**The primary source of k_i, k_a, k_s, k_r_max and K_D**] <https://pubmed.ncbi.nlm.nih.gov/15498514/>
10. Worek F, Eyer P, Aurbek N, Szinicz L, Thiermann H. Recent advances in evaluation of oxime efficacy in nerve agent poisoning by in vitro analysis. *Toxicol Appl Pharmacol* 2007;219:226-34. [Formalisation of the saturating reactivation expression] <https://pubmed.ncbi.nlm.nih.gov/17169391/>
11. Aurbek N, Thiermann H, Szinicz L, Eyer P, Worek F. Analysis of inhibition, reactivation and aging kinetics of highly toxic organophosphorus compounds with human and pig acetylcholinesterase. *Toxicology* 2006;224:91-9. [Ageing half-life diethyl 33 h / dimethyl 3.7 h] <https://pubmed.ncbi.nlm.nih.gov/16720069/>
12. Worek F, Diepold C, Eyer P. Dimethylphosphoryl-inhibited human cholinesterases: inhibition, reactivation, and aging kinetics. *Arch Toxicol* 1999;73:7-14. [The fast spontaneous reactivation of dimethyl OPs, t½ 0.7 h — the reason the model computes φ *lower* than for diethyl] <https://pubmed.ncbi.nlm.nih.gov/10207609/>
13. Worek F, Thiermann H, Wille T. Oximes in organophosphate poisoning: 60 years of hope and despair. *Chem Biol Interact* 2016;259:93-8. [A mechanistic account of why oximes fail] <https://pubmed.ncbi.nlm.nih.gov/27125983/>
14. Eyer P. The role of oximes in the management of organophosphorus pesticide poisoning. *Toxicol Rev* 2003;22:165-90. [Review of oxime pharmacokinetics and pharmacodynamics · the re-inhibition phenomenon] <https://pubmed.ncbi.nlm.nih.gov/15181665/>
15. Mason HJ. The recovery of plasma cholinesterase and erythrocyte acetylcholinesterase activity in workers after over-exposure to dichlorvos. *Occup Med (Lond)* 2000;50:343-7. [Erythrocyte AChE is not resynthesised — the `KRBCNEW` value] <https://pubmed.ncbi.nlm.nih.gov/10975133/>
16. Thiermann H, Szinicz L, Eyer P, et al. Correlation between red blood cell acetylcholinesterase activity and neuromuscular transmission in organophosphate poisoning. *Chem Biol Interact* 2005;157-158:345-7. [The correlation between RBC AChE and neuromuscular transmission — and its limits] <https://pubmed.ncbi.nlm.nih.gov/16289002/>
17. Worek F, Mast U, Kiderlen D, Diepold C, Eyer P. Improved determination of acetylcholinesterase activity in human whole blood. *Clin Chim Acta* 1999;288:73-90. [The AChE assay — the definition of the biomarker outputs] <https://pubmed.ncbi.nlm.nih.gov/10529460/>

---

## 3. Toxicokinetics and bioactivation

In the model, **what inhibits AChE is not the parent compound but the oxon that CYP makes**, and
the fact that its formation saturates (`VMAXBIO`, `KMBIO`) is what produces the conclusion that "at
large ingested doses the oxon concentration becomes independent of dose".

18. Eyer F, Meischner V, Kiderlen D, et al. Human parathion poisoning: a toxicokinetic analysis. *Toxicol Rev* 2003;22:143-63. [Measured blood concentrations of parent compound and oxon; redistribution from the fat depot] <https://pubmed.ncbi.nlm.nih.gov/15181664/>
19. Eyer F, Roberts DM, Buckley NA, et al. Extreme variability in the formation of chlorpyrifos oxon (CPO) in patients poisoned by chlorpyrifos (CPF). *Biochem Pharmacol* 2009;78:531-7. [**The extreme between-patient variability of oxon formation — the basis for `ETA_BIO` and `PON1SC` in the model, and the most uncertain input to Ω**] <https://pubmed.ncbi.nlm.nih.gov/19433067/>
20. Timchalk C, Nolan RJ, Mendrala AL, et al. A physiologically based pharmacokinetic and pharmacodynamic (PBPK/PD) model for the organophosphate insecticide chlorpyrifos in rats and humans. *Toxicol Sci* 2002;66:34-53. [The chlorpyrifos PBPK structure · the fat compartment] <https://pubmed.ncbi.nlm.nih.gov/11861971/>
21. Buratti FM, Volpe MT, Meneguz A, Vittozzi L, Testai E. CYP-specific bioactivation of four organophosphorothioate pesticides by human liver microsomes. *Toxicol Appl Pharmacol* 2003;186:143-54. [The CYP2B6/3A4/1A2 branch point between desulfuration and dearylation] <https://pubmed.ncbi.nlm.nih.gov/12620367/>
22. Costa LG, Cole TB, Vitalone A, Furlong CE. Measurement of paraoxonase (PON1) status as a potential biomarker of susceptibility to organophosphate toxicity. *Clin Chim Acta* 2005;352:37-47. [PON1 Q192R — the basis for `PON1SC = 0.33`] <https://pubmed.ncbi.nlm.nih.gov/15653098/>
23. Furlong CE, Marsillach J, Jarvik GP, Costa LG. Paraoxonases-1, -2 and -3: what are their functions? *Chem Biol Interact* 2016;259:51-62. [PON1 substrate specificity: dimethyl oxons are barely hydrolysed] <https://pubmed.ncbi.nlm.nih.gov/27238723/>
24. Eddleston M, Worek F, Eyer P, et al. Poisoning with the S-alkyl organophosphorus insecticides profenofos and prothiofos. *QJM* 2009;102:785-92. [The kinetics of atypical OPs — the basis for extending the model's library] <https://pubmed.ncbi.nlm.nih.gov/19737791/>
25. Eddleston M, Street JM, Self I, et al. A role for solvents in the toxicity of agricultural organophosphorus pesticides. *Toxicology* 2012;294:94-103. [**The independent toxicity of the solvents (cyclohexanone · hydrocarbons) — the model's `SOLVCV` and `SOLVRESP`**] <https://pubmed.ncbi.nlm.nih.gov/22365945/>
26. Roberts DM, Aaron CK. Management of acute organophosphorus pesticide poisoning. *BMJ* 2007;334:629-34. [The absorption window, and the time dependence of activated charcoal] <https://pubmed.ncbi.nlm.nih.gov/17379909/>

---

## 4. Atropine — a closed-loop drug

In the model atropine is implemented **as a controller rather than as a dosing schedule**, with the
consequence that the cumulative dose becomes an *output* of the model. The references below are the
basis for its targets (a dry chest, HR > 80, systolic blood pressure > 80) and for the doubling protocol.

27. Eddleston M, Buckley NA, Checketts H, et al. Speed of initial atropinisation in significant organophosphorus pesticide poisoning — a systematic comparison of recommended regimens. *J Toxicol Clin Toxicol* 2004;42:865-75. [The doubling protocol vs conventional titration — the model's `ATRMODE`] <https://pubmed.ncbi.nlm.nih.gov/15533024/>
28. Eddleston M, Dawson A, Karalliedde L, et al. Early management after self-poisoning with an organophosphorus or carbamate pesticide — a treatment protocol for junior doctors. *Crit Care* 2004;8:R391-7. [An operational definition of the atropinisation endpoint] <https://pubmed.ncbi.nlm.nih.gov/15566582/>
29. Abedin MJ, Sayeed AA, Basher A, et al. Open-label randomized clinical trial of atropine bolus injection versus incremental boluses plus infusion for organophosphate poisoning in Bangladesh. *J Med Toxicol* 2012;8:108-17. [Infusion vs bolus — `KP_RAPID`/`KP_SLOW`] <https://pubmed.ncbi.nlm.nih.gov/22331327/>
30. Perera PMS, Shahmy S, Gawarammana I, Dawson AH. Comparison of two commonly practiced atropinization regimens in acute organophosphorus and carbamate poisoning. *Hum Exp Toxicol* 2008;27:513-8. [The measured range of cumulative atropine dose] <https://pubmed.ncbi.nlm.nih.gov/18784205/>
31. Ali-Melkkilä T, Kanto J, Iisalo E. Pharmacokinetics and related pharmacodynamics of anticholinergic drugs. *Acta Anaesthesiol Scand* 1993;37:633-42. [**Atropine V_d and CL, and the quaternary ammonium structure of glycopyrrolate → no BBB penetration: the model's `BBBF` 0.35 vs 0.02**] <https://pubmed.ncbi.nlm.nih.gov/8249551/>
32. Bardin PG, van Eeden SF. Organophosphate poisoning: grading the severity and comparing treatment between atropine and glycopyrrolate. *Crit Care Med* 1990;18:956-60. [The glycopyrrolate comparison — scenario S15] <https://pubmed.ncbi.nlm.nih.gov/2394120/>
33. Connors NJ, Harnett ZH, Hoffman RS. Comparison of current recommended regimens of atropinization in organophosphate poisoning. *J Med Toxicol* 2014;10:143-7. [A comparison between the recommendations] <https://pubmed.ncbi.nlm.nih.gov/24619543/>

---

## 5. Oximes — pharmacokinetics and the trials

34. Eddleston M, Eyer P, Worek F, et al. Pralidoxime in acute organophosphorus insecticide poisoning — a randomised controlled trial. *PLoS Med* 2009;6:e1000104. [**The randomised trial of the WHO regimen (30 mg/kg + 8 mg/kg/h), with a negative result — what the virtual trial in §13 of the model sets out to reproduce**] <https://pubmed.ncbi.nlm.nih.gov/19564902/>
35. Buckley NA, Eddleston M, Li Y, Bevan M, Robertson J. Oximes for acute organophosphate pesticide poisoning. *Cochrane Database Syst Rev* 2011;(2):CD005085. [Systematic review: no benefit] <https://pubmed.ncbi.nlm.nih.gov/21328273/>
36. Pawar KS, Bhoite RR, Pillay CP, Chavan SC, Malshikare DS, Garad SG. Continuous pralidoxime infusion versus repeated bolus injection to treat organophosphorus pesticide poisoning: a randomised controlled trial. *Lancet* 2006;368:2136-41. [The superiority of continuous infusion — the model's S02 vs S03 contrast] <https://pubmed.ncbi.nlm.nih.gov/17174705/>
37. Medicis JJ, Stork CM, Howland MA, Hoffman RS, Goldfrank LR. Pharmacokinetics following a loading plus a continuous infusion of pralidoxime compared with the traditional short infusion regimen in human volunteers. *J Toxicol Clin Toxicol* 1996;34:289-95. [The pralidoxime PK parameters `OXV1/OXCL`] <https://pubmed.ncbi.nlm.nih.gov/8667465/>
38. Thiermann H, Szinicz L, Eyer F, et al. Modern strategies in therapy of organophosphate poisoning. *Toxicol Lett* 1999;107:233-9. [The obidoxime regimen 250 mg + 750 mg/24 h] <https://pubmed.ncbi.nlm.nih.gov/10414801/>
39. Thiermann H, Zilker T, Eyer F, Felgenhauer N, Eyer P, Worek F. Monitoring of neuromuscular transmission in organophosphate pesticide-poisoned patients. *Toxicol Lett* 2009;191:297-304. [The conditions under which reactivation is observed in real patients] <https://pubmed.ncbi.nlm.nih.gov/19825404/>
40. Thiermann H, Worek F, Kehe K. Limitations and challenges in treatment of acute chemical warfare agent poisoning. *Chem Biol Interact* 2013;206:435-43. [An earlier qualitative statement of the oxime-ceiling concept] <https://pubmed.ncbi.nlm.nih.gov/23994499/>
41. Kiderlen D, Eyer P, Worek F. Formation and disposition of diethylphosphoryl-obidoxime, a potent anticholinesterase that is hydrolyzed by human paraoxonase (PON1). *Biochem Pharmacol* 2005;69:1853-67. [Re-inhibition by the phosphorylated oxime — the model's `REBOUND` node] <https://pubmed.ncbi.nlm.nih.gov/15869746/>
42. Rahimi R, Nikfar S, Abdollahi M. Increased morbidity and mortality in acute human organophosphate-poisoned patients treated by oximes: a meta-analysis of clinical trials. *Hum Exp Toxicol* 2006;25:157-62. [**The meta-analysis in which mortality was higher in the oxime-treated group — consistent with the direction the model (unintentionally) reproduces, and flagged as the 'most exposed prediction' in the README**] <https://pubmed.ncbi.nlm.nih.gov/16696287/>
43. Peter JV, Moran JL, Graham P. Oxime therapy and outcomes in human organophosphate poisoning: an evaluation using meta-analytic techniques. *Crit Care Med* 2006;34:502-10. [An independent meta-analysis on the same question] <https://pubmed.ncbi.nlm.nih.gov/16424734/>
44. Syed S, Gurcoo SA, Farooqui AK, Nisa W, Sofi K, Wani TM. Is the World Health Organization-recommended dose of pralidoxime effective in the treatment of organophosphorus poisoning? *Indian J Anaesth* 2015;59:31-6. [An attempt to re-test the WHO regimen] <https://pubmed.ncbi.nlm.nih.gov/25684811/>

---

## 6. The nicotinic limb and the intermediate syndrome

45. Senanayake N, Karalliedde L. Neurotoxic effects of organophosphorus insecticides: an intermediate syndrome. *N Engl J Med* 1987;316:761-3. [The first description of the intermediate syndrome — the model's `Rn_des`/`IMS`] <https://pubmed.ncbi.nlm.nih.gov/3029588/>
46. Jayawardane P, Dawson AH, Weerasinghe V, Karalliedde L, Buckley NA, Senanayake N. The spectrum of intermediate syndrome following acute organophosphate poisoning: a prospective cohort study from Sri Lanka. *PLoS Med* 2008;5:e147. [Onset at 24-96 h, duration 4-18 days — `KDES`, `KRECDES`] <https://pubmed.ncbi.nlm.nih.gov/18630983/>
47. Jayawardane P, Senanayake N, Buckley NA, Dawson AH. Electrophysiological correlates of the intermediate syndrome in acute organophosphate poisoning. *Clin Toxicol* 2012;50:250-3. [Decrement on repetitive stimulation — the electrophysiology of depolarising blockade] <https://pubmed.ncbi.nlm.nih.gov/22385107/>
48. Karalliedde L, Baker D, Marrs TC. Organophosphate-induced intermediate syndrome: aetiology and relationships with myopathy. *Toxicol Rev* 2006;25:1-14. [A survey of the mechanistic hypotheses] <https://pubmed.ncbi.nlm.nih.gov/16856766/>
49. Wood SJ, Slater CR. Safety factor at the neuromuscular junction. *Prog Neurobiol* 2001;64:393-429. [**The safety factor at the neuromuscular junction — the model's `SF_NMJ` = 0.28, i.e. the structure in which resting ventilation requires only 28% of maximal inspiratory force**] <https://pubmed.ncbi.nlm.nih.gov/11275359/>
50. Yang CC, Deng JF. Intermediate syndrome following organophosphate insecticide poisoning. *J Chin Med Assoc* 2007;70:467-72. [The clinical course and the requirement for mechanical ventilation] <https://pubmed.ncbi.nlm.nih.gov/18063499/>

---

## 7. Respiratory failure and supportive care

More than fifty clinical observations point one way: **the single factor that decides survival in
this disease is the availability of mechanical ventilation**, and scenario S13 of the model (no
ventilator) quantifies it.

51. Eddleston M, Mohamed F, Davies JOJ, et al. Respiratory failure in acute organophosphorus pesticide self-poisoning. *QJM* 2006;99:513-22. [**The decomposition of respiratory failure into central depression · secretions · neuromuscular blockade — the direct basis for the structure of §8 of the model**] <https://pubmed.ncbi.nlm.nih.gov/16829539/>
52. Hulse EJ, Davies JOJ, Simpson AJ, Sciuto AM, Eddleston M. Respiratory complications of organophosphorus nerve agent and insecticide poisoning: implications for respiratory and critical care. *Am J Respir Crit Care Med* 2014;190:1342-54. [Airway secretions, aspiration, ARDS] <https://pubmed.ncbi.nlm.nih.gov/25419614/>
53. Gaspari RJ, Paydarfar D. Respiratory recovery following organophosphate poisoning in a rat model is suppressed by isolated hypoxia at the point of apnea. *Toxicology* 2011;290:59-64. [The hypoxic vicious circle of central apnoea] <https://pubmed.ncbi.nlm.nih.gov/21875643/>
54. Carey JL, Dunn C, Gaspari RJ. Central respiratory failure during acute organophosphate poisoning. *Respir Physiol Neurobiol* 2013;189:403-10. [Failure of the preBötzinger rhythm generator — the model's `PREBOT`/`DMAXNONR`] <https://pubmed.ncbi.nlm.nih.gov/23933009/>
55. Eddleston M, Chowdhury FR. Pharmacological treatment of organophosphorus insecticide poisoning: the old and the (possible) new. *Br J Clin Pharmacol* 2016;81:462-70. [A critical survey of current and candidate treatments] <https://pubmed.ncbi.nlm.nih.gov/26366467/>
56. Munidasa UADD, Gawarammana IB, Kularatne SAM, Kumarasiri PVR, Goonasekera CDA. Survival pattern in patients with acute organophosphate poisoning receiving intensive care. *J Toxicol Clin Toxicol* 2004;42:343-7. [Intensive-care survival curves — the calibration standard for the model's mortality] <https://pubmed.ncbi.nlm.nih.gov/15461240/>
57. Senanayake N, de Silva HJ, Karalliedde L. A scale to assess severity in organophosphorus intoxication: POP scale. *Hum Exp Toxicol* 1993;12:297-9. [**The Peradeniya OP Poisoning scale — the score reconstructed from the state vector in the model's `$TABLE`**] <https://pubmed.ncbi.nlm.nih.gov/8104007/>
58. Davies JOJ, Eddleston M, Buckley NA. Predicting outcome in acute organophosphorus poisoning with a poison severity score or the Glasgow coma scale. *QJM* 2008;101:371-9. [The predictive power of the severity score] <https://pubmed.ncbi.nlm.nih.gov/18319295/>

---

## 8. Central nervous system · seizures · benzodiazepines

59. McDonough JH Jr, Shih TM. Neuropharmacological mechanisms of nerve agent-induced seizure and neuropathology. *Neurosci Biobehav Rev* 1997;21:559-79. [The cholinergic → glutamatergic transition, and the 30-minute window — the model's `SEIZ`→`NMDA`] <https://pubmed.ncbi.nlm.nih.gov/9353792/>
60. Shih TM, Duniho SM, McDonough JH. Control of nerve agent-induced seizures is critical for neuroprotection and survival. *Toxicol Appl Pharmacol* 2003;188:69-80. [The time-dependent effect of diazepam] <https://pubmed.ncbi.nlm.nih.gov/12691725/>
61. Eddleston M, Buckley NA, Checketts H, et al. Diazepam in organophosphorus poisoning. (In: Management reviews above; see also) Marrs TC, Rice P, Vale JA. The role of oximes in the treatment of nerve agent poisoning in civilian casualties. *Toxicol Rev* 2006;25:297-323. [The basis for co-administering diazepam] <https://pubmed.ncbi.nlm.nih.gov/17288500/>
62. Bird SB, Gaspari RJ, Dickson EW. Early death due to severe organophosphate poisoning is a centrally mediated process. *Acad Emerg Med* 2003;10:295-8. [**Early death is centrally mediated — the reason the model connects the BBB penetration of atropine to outcome**] <https://pubmed.ncbi.nlm.nih.gov/12670843/>
63. Chen Y. Organophosphate-induced brain damage: mechanisms, neuropsychiatric and neurological consequences, and potential therapeutic strategies. *Neurotoxicology* 2012;33:391-400. [COPIND] <https://pubmed.ncbi.nlm.nih.gov/22498093/>

---

## 9. Delayed polyneuropathy (OPIDN) and NTE

64. Lotti M, Moretto A. Organophosphate-induced delayed polyneuropathy. *Toxicol Rev* 2005;24:37-49. [**NTE inhibition plus ageing above 70% is the condition for onset — the model's `NTETH`**] <https://pubmed.ncbi.nlm.nih.gov/16042503/>
65. Glynn P. Neuropathy target esterase and phospholipid deacylation. *Biochim Biophys Acta* 2005;1736:87-93. [The biochemistry of NTE (PNPLA6)] <https://pubmed.ncbi.nlm.nih.gov/16137924/>
66. Jokanović M, Kosanović M, Brkić D, Vukomanović P. Organophosphate induced delayed polyneuropathy in man: an overview. *Clin Neurol Neurosurg* 2011;113:7-10. [The list of causative compounds — dimethyl OPs barely cause it] <https://pubmed.ncbi.nlm.nih.gov/20880635/>

---

## 10. Adjunctive treatments · the failed candidates

67. Pajoumand A, Shadnia S, Rezaie A, Abdi M, Abdollahi M. Benefits of magnesium sulfate in the management of acute human poisoning by organophosphorus insecticides. *Hum Exp Toxicol* 2004;23:565-9. [Magnesium sulfate — scenario S16] <https://pubmed.ncbi.nlm.nih.gov/15688984/>
68. Basher A, Rahman SH, Ghose A, Arif SM, Faiz MA, Dawson AH. Phase II study of magnesium sulfate in acute organophosphate pesticide poisoning. *Clin Toxicol* 2013;51:35-40. [Dose finding] <https://pubmed.ncbi.nlm.nih.gov/23148565/>
69. Perera PM, Jayamanna SF, Hettiarachchi R, et al. A phase II clinical trial to assess the safety of clonidine in acute organophosphorus pesticide poisoning. *Trials* 2009;10:73. [Clonidine — the model's `CLONID` node] <https://pubmed.ncbi.nlm.nih.gov/19678923/>
70. Eddleston M, Juszczak E, Buckley NA, et al. Multiple-dose activated charcoal in acute self-poisoning: a randomised controlled trial. *Lancet* 2008;371:579-87. [**Multiple-dose activated charcoal: no benefit — the model offers the arithmetic of the absorption window as the reason**] <https://pubmed.ncbi.nlm.nih.gov/18280328/>
71. Li Y, Tse ML, Gawarammana I, Buckley N, Eddleston M. Systematic review of controlled clinical trials of gastric lavage in acute organophosphorus pesticide poisoning. *Clin Toxicol* 2009;47:179-92. [Gastric lavage] <https://pubmed.ncbi.nlm.nih.gov/19306191/>
72. Peter JV, Moran JL, Graham PL. Advances in the management of organophosphate poisoning. *Expert Opin Pharmacother* 2007;8:1451-64. [A review of adjunctive therapy] <https://pubmed.ncbi.nlm.nih.gov/17661728/>
73. Chowdhary S, Bhattacharyya R, Banerjee D. Acute organophosphorus poisoning. *Clin Chim Acta* 2014;431:66-76. [The laboratory-medicine perspective] <https://pubmed.ncbi.nlm.nih.gov/24508989/>

---

## 11. Biocatalysts / bioscavengers — why they do not work for insecticides

Section 6 of the model shows that **a stoichiometric scavenger is arithmetically impossible in
insecticide self-poisoning** (mmol ingested versus the µmol of scavenger that can be given). The same
strategy does hold for nerve-agent prophylaxis, and that the difference lies precisely in the molar amount is consistent with the references below.

74. Lenz DE, Yeung D, Smith JR, Sweeney RE, Lumley LA, Cerasoli DM. Stoichiometric and catalytic scavengers as protection against nerve agent toxicity: a mini review. *Toxicology* 2007;233:31-9. [The principle of stoichiometric scavenging and the molar constraint] <https://pubmed.ncbi.nlm.nih.gov/17188793/>
75. Nachon F, Brazzolotto X, Trovaslet M, Masson P. Progress in the development of enzyme-based nerve agent bioscavengers. *Chem Biol Interact* 2013;206:536-44. [Catalytic scavengers] <https://pubmed.ncbi.nlm.nih.gov/23811386/>
76. Rosenberg YJ, Laube B, Mao L, et al. Pulmonary delivery of an aerosolized recombinant human butyrylcholinesterase pretreatment protects against aerosolized paraoxon in macaques. *Chem Biol Interact* 2013;203:167-71. [rHuBChE pretreatment] <https://pubmed.ncbi.nlm.nih.gov/23103586/>
77. Ashani Y, Pistinner S. Estimation of the upper limit of human butyrylcholinesterase dose required for protection against organophosphates toxicity: a mathematically based toxicokinetic model. *Toxicol Sci* 2004;77:358-67. [**The mathematical upper bound on the required BChE dose — the calculation that anticipates §6 of the model**] <https://pubmed.ncbi.nlm.nih.gov/14691206/>

---

## 12. Cardiovascular toxicity

78. Ludomirsky A, Klein HO, Sarelli P, et al. Q-T prolongation and polymorphous ("torsade de pointes") ventricular arrhythmias associated with organophosphorus insecticide poisoning. *Am J Cardiol* 1982;49:1654-8. [QTc · torsade de pointes] <https://pubmed.ncbi.nlm.nih.gov/7081053/>
79. Yurumez Y, Yavuz Y, Saglam H, et al. Electrocardiographic findings of acute organophosphate poisoning. *J Emerg Med* 2009;36:39-42. [The frequency of the ECG findings] <https://pubmed.ncbi.nlm.nih.gov/18024061/>
80. Chuang FR, Jang SW, Lin JL, Chern MS, Chen JB, Hsu KT. QTc prolongation indicates a poor prognosis in patients with organophosphate poisoning. *Am J Emerg Med* 1996;14:451-3. [A prognostic marker] <https://pubmed.ncbi.nlm.nih.gov/8765108/>
81. Anand S, Singh S, Nahar Saikia U, Bhalla A, Paul Sharma Y, Singh D. Cardiac abnormalities in acute organophosphate poisoning. *Clin Toxicol* 2009;47:230-5. [Myocardial injury] <https://pubmed.ncbi.nlm.nih.gov/19306191/>

---

## 13. Modelling methodology (QSP / PBPK)

82. Baron KT, Gastonguay MR. Simulation from ODE-based population PK/PD and systems pharmacology models in R with mrgsolve. *J Pharmacokinet Pharmacodyn* 2015;42:S84-5. [mrgsolve] <https://mrgsolve.org>
83. Timchalk C, Poet TS. Development of a physiologically based pharmacokinetic and pharmacodynamic model to determine dosimetry and cholinesterase inhibition for a binary mixture of chlorpyrifos and diazinon in the rat. *Neurotoxicology* 2008;29:428-43. [A precursor OP PBPK/PD model] <https://pubmed.ncbi.nlm.nih.gov/18359092/>
84. Gearhart JM, Jepson GW, Clewell HJ 3rd, Andersen ME, Conolly RB. Physiologically based pharmacokinetic and pharmacodynamic model for the inhibition of acetylcholinesterase by diisopropylfluorophosphate. *Toxicol Appl Pharmacol* 1990;106:295-310. [The prototype PBPK/PD model of AChE inhibition] <https://pubmed.ncbi.nlm.nih.gov/2260096/>
85. Bhattacharya S, Conolly RB, Kaminski NE, Thomas RS, Andersen ME, Zhang Q. A bistable switch underlying B-cell differentiation and its disruption by the environmental contaminant TCDD. *Toxicol Sci* 2010;115:51-65. [A methodology for analysing switches in toxicological systems — the conceptual frame used to analyse this model's three-state switch] <https://pubmed.ncbi.nlm.nih.gov/20123757/>
86. Sterner TR, Ruark CD, Covington TR, Yu KO, Gearhart JM. A physiologically based pharmacokinetic model for the oxime TMB-4: simulation of rodent and human data. *Arch Toxicol* 2013;87:661-80. [Oxime PBPK] <https://pubmed.ncbi.nlm.nih.gov/23151815/>

---

## 14. Guidelines and general reviews

87. World Health Organization. *The WHO Recommended Classification of Pesticides by Hazard and Guidelines to Classification 2019.* Geneva: WHO, 2020. [The WHO Class I/II classification — the model's `REG` node] <https://www.who.int/publications/i/item/9789240005662>
88. Eddleston M. Novel clinical toxicology and pharmacology of organophosphorus insecticide self-poisoning. *Annu Rev Pharmacol Toxicol* 2019;59:341-60. [**The best single summary of the whole clinical narrative this model sets out to reproduce**] <https://pubmed.ncbi.nlm.nih.gov/30110584/>
89. Vale A, Lotti M. Organophosphorus and carbamate insecticide poisoning. *Handb Clin Neurol* 2015;131:149-68. [A review from the neurological perspective] <https://pubmed.ncbi.nlm.nih.gov/26563788/>
90. Robb EL, Baker MB. Organophosphate toxicity. *StatPearls* [Internet]. Treasure Island (FL): StatPearls Publishing, 2023. [An educational summary] <https://www.ncbi.nlm.nih.gov/books/NBK470430/>

---

## Appendix A — Provenance table for each model parameter

| Parameter | Value (default: chlorpyrifos) | Source |
|---|---|---|
| `KI` AChE inhibition constant of the oxon | 0.18 nM⁻¹h⁻¹ (= 3.0×10⁶ M⁻¹min⁻¹) | 9, 11 |
| `T12AGE` ageing half-life | diethyl 33 h / dimethyl 3.7 h | 9, 11, 12 |
| `T12SPO` spontaneous reactivation half-life | diethyl 77 h / dimethyl 0.7 h | 9, 12 |
| `KRMAX`, `KDOX` saturation of oxime reactivation | 2-PAM 36 h⁻¹ / 200 µM; obidoxime 48 h⁻¹ / 40 µM | 9, 10, 13 |
| `VMAXBIO`, `KMBIO` CYP desulfuration | 6.0×10⁴ nmol/h, 20 µM | 19, 20, 21 |
| `CLOXON`, `PON1SC` oxon hydrolysis | 300 L/h, ×0.33 in QQ | 22, 23 |
| `VFTH` fat depot | 700 L (fenthion 1750 L) | 18, 20 |
| `OXV1/OXCL` pralidoxime PK | 24 L / 21 L/h | 34, 37 |
| `ATR_V/ATR_CL/BBBF` atropine | 200 L / 46 L/h / 0.35 | 31 |
| `BBBF` glycopyrrolate | 0.02 | 31, 32 |
| `SF_NMJ` neuromuscular safety factor | 0.28 | 49 |
| `KDES/KRECDES` intermediate syndrome | 0.075 / 0.050 h⁻¹ | 45, 46, 47 |
| `NTETH` OPIDN threshold | 70% aged NTE | 64, 66 |
| `HVENT` risk of a ventilation complication | 0.0012 h⁻¹ (≈2.9%/day) | 51, 52, 56 |
| `SOLVCV` myocardial depression by the solvent | cyclohexanone 1.8 vs hydrocarbons 0.4-0.5 | 25 |
| Distribution of ingested volume in the virtual trial | log-normal, median 28 mL | 4, 6, 34 |

## Appendix B — Clinical observations the model does and does not reproduce

| Clinical observation | Model | Supporting reference |
|---|---|---|
| Dimethyl OPs have a higher case fatality than diethyl | Reproduced (ageing rate + oxon concentration) | 4 |
| Oximes reactivate AChE but do not change mortality | Reproduced (§13 virtual trial RR ≈ 1) | 34, 35, 43 |
| Continuous infusion is better than bolus | Reproduced (S02 vs S03) | 36 |
| Rapid atropine doubling is better than conventional titration | Reproduced (S02 vs S14) | 27, 29 |
| Multiple-dose activated charcoal gives no benefit | Reproduced (S17 vs S18: the effect is lost beyond one hour) | 70 |
| Early death is centrally mediated | Reproduced (worse with glycopyrrolate in S15) | 62 |
| BChE is a prognostic marker | Partly reproduced (provided as an output; the prognostic equation is not fitted) | 5 |
| The meta-analysis in which mortality is in fact higher in the oxime group | **Reproduced unintentionally** — see the 'most exposed prediction' in the README | 42, 43 |
