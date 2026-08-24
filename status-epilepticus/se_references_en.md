# Status epilepticus QSP model — references
# Status Epilepticus QSP Model — References

Every structural claim and parameter of this model corresponds to a reference below.
**Every PMID was looked up with NCBI E-utilities and its author, journal, year, and title confirmed.**

The `→` mark at the end of each entry indicates **which part** of the model that reference was used in.
(`§3` = the cluster number in the mechanistic map, `PARAM` = an mrgsolve parameter, `CAL` = a calibration basis)

---

## 1. Definition · epidemiology · guidelines

1. Trinka E, et al. **A definition and classification of status epilepticus — Report of the ILAE Task Force on Classification of Status Epilepticus.** Epilepsia 2015. — [PMID 26336950](https://pubmed.ncbi.nlm.nih.gov/26336950/)
   → t1 = 5 minutes (the transition to self-sustainment), t2 = 30 minutes (the onset of injury). The time scale of the model's `N_transit` node and `INJURY` integrator comes from here. §2, §9

2. Glauser T, et al. **Evidence-Based Guideline: Treatment of Convulsive Status Epilepticus in Children and Adults (American Epilepsy Society).** Epilepsy Curr 2016. — [PMID 26900382](https://pubmed.ncbi.nlm.nih.gov/26900382/)
   → the reference values for the timing and dose of the first, second, and third line. Scenarios 02 and 21

3. Brophy GM, et al. **Guidelines for the evaluation and management of status epilepticus.** Neurocrit Care 2012. — [PMID 22528274](https://pubmed.ncbi.nlm.nih.gov/22528274/)
   → the dose range of the third-line anaesthetics (midazolam 0.05-2 mg/kg/h, propofol 2-10 mg/kg/h). §14, PARAM

4. DeLorenzo RJ, et al. **A prospective, population-based epidemiologic study of status epilepticus in Richmond, Virginia.** Neurology 1996. — [PMID 8780085](https://pubmed.ncbi.nlm.nih.gov/8780085/)
   → incidence and the distribution by aetiology. §1

5. Betjemann JP, Lowenstein DH. **Status epilepticus in adults.** Lancet Neurol 2015. — [PMID 25908090](https://pubmed.ncbi.nlm.nih.gov/25908090/)
   → the overall clinical frame. All clusters

6. Neligan A, Shorvon SD. **Frequency and prognosis of convulsive status epilepticus of different causes: a systematic review.** Arch Neurol 2010. — [PMID 20697043](https://pubmed.ncbi.nlm.nih.gov/20697043/)
   → mortality by aetiology — the basis for dissociation 3 of §18, that `EDRIVE` dominates the variance in outcome

7. Neligan A, et al. **Change in mortality of generalized convulsive status epilepticus in high-income countries over time.** JAMA Neurol 2019. — [PMID 31135807](https://pubmed.ncbi.nlm.nih.gov/31135807/)
   → the observation that mortality has not changed much despite advances in treatment — consistent with the model's claim that "the clock dominates the choice"

8. Kantanen AM, et al. **Incidence and mortality of super-refractory status epilepticus in adults.** Epilepsy Behav 2015. — [PMID 26141934](https://pubmed.ncbi.nlm.nih.gov/26141934/)
   → the frequency and mortality of SRSE. §18 `K_rse`, `K_mort`

9. Kantanen AM, et al. **Long-term outcome of refractory status epilepticus in adults: a retrospective population-based study.** Epilepsy Res 2017. — [PMID 28402834](https://pubmed.ncbi.nlm.nih.gov/28402834/)

10. Kantanen AM, et al. **Predictors of hospital and one-year mortality in intensive care patients with refractory status epilepticus.** Crit Care 2017. — [PMID 28330483](https://pubmed.ncbi.nlm.nih.gov/28330483/)

11. Lowenstein DH, Bleck T, Macdonald RL. **It's time to revise the definition of status epilepticus.** Epilepsia 1999. — [PMID 9924914](https://pubmed.ncbi.nlm.nih.gov/9924914/)
    → the argument for the shift from a 30-minute to a 5-minute definition is exactly the clock logic of this model

12. Fountain NB. **Status epilepticus: risk factors and complications.** Epilepsia 2000. — [PMID 10885737](https://pubmed.ncbi.nlm.nih.gov/10885737/) → §10, §11

13. Mayer SA, et al. **Refractory status epilepticus: frequency, risk factors, and impact on outcome.** Arch Neurol 2002. — [PMID 11843690](https://pubmed.ncbi.nlm.nih.gov/11843690/)

14. Novy J, Logroscino G, Rossetti AO. **Refractory status epilepticus: a prospective observational study.** Epilepsia 2010. — [PMID 19817823](https://pubmed.ncbi.nlm.nih.gov/19817823/)

---

## 2. Self-sustainment and time dependence

15. Lothman E. **The biochemical basis and pathophysiology of status epilepticus.** Neurology 1990. — [PMID 2185436](https://pubmed.ncbi.nlm.nih.gov/2185436/)
    → the transition from the compensated to the decompensated phase (about 30 minutes). §10 → §11

16. Mazarati AM, et al. **Self-sustaining status epilepticus after brief electrical stimulation of the perforant path.** Brain Res 1998. — [PMID 9729413](https://pubmed.ncbi.nlm.nih.gov/9729413/)
    → the seizure continues even when the provoking stimulus is removed = the bistability of the `SEIZ` state variable. §2

17. Wasterlain CG, et al. **Self-sustaining status epilepticus: a condition maintained by potentiation of glutamate receptors and by plasticity of the GABA-A receptor.** Epilepsia 2000. — [PMID 10999535](https://pubmed.ncbi.nlm.nih.gov/10999535/)
    → **the source of the two-clock concept of this model itself.** §3, §5

18. Wasterlain CG, Chen JW. **Mechanistic and pharmacologic aspects of status epilepticus and its treatment with new antiepileptic drugs.** Epilepsia 2008. — [PMID 19087119](https://pubmed.ncbi.nlm.nih.gov/19087119/)

19. Chen JW, Wasterlain CG. **Status epilepticus: pathophysiology and management in adults.** Lancet Neurol 2006. — [PMID 16488380](https://pubmed.ncbi.nlm.nih.gov/16488380/)

20. Lowenstein DH, Alldredge BK. **Status epilepticus.** N Engl J Med 1998. — [PMID 9521986](https://pubmed.ncbi.nlm.nih.gov/9521986/)

21. Shorvon S, Ferlisi M. **The treatment of super-refractory status epilepticus: a critical review of available therapies and a clinical treatment protocol.** Brain 2011. — [PMID 21914716](https://pubmed.ncbi.nlm.nih.gov/21914716/) → §14

22. Ferlisi M, Shorvon S. **The outcome of therapies in refractory and super-refractory convulsive status epilepticus and recommendations for therapy.** Brain 2012. — [PMID 22577217](https://pubmed.ncbi.nlm.nih.gov/22577217/)

23. Bleck TP. **Refractory status epilepticus.** Curr Opin Crit Care 2005. — [PMID 15758590](https://pubmed.ncbi.nlm.nih.gov/15758590/)

---

## 3. **Clock 1** — internalisation of the synaptic γ2-GABA-A receptor (the central axis of the model)

24. **Kapur J, Macdonald RL. Rapid seizure-induced reduction of benzodiazepine and Zn²⁺ sensitivity of hippocampal dentate granule cell GABA-A receptors.** J Neurosci 1997. — [PMID 9295398](https://pubmed.ncbi.nlm.nih.gov/9295398/)
    → **the potency of diazepam falls about 20-fold against about 3-fold for phenobarbital.** The model does not take these numbers as inputs; it reproduces them from the three pools (§3, §4) as `CREQ_BZD` / `CREQ_PB` (7.0-fold vs 2.6-fold @30 min — direction and asymmetry reproduced, magnitude underestimated). CAL

25. **Naylor DE, Liu H, Wasterlain CG. Trafficking of GABA(A) receptors, loss of inhibition, and a mechanism for pharmacoresistance in status epilepticus.** J Neurosci 2005. — [PMID 16120773](https://pubmed.ncbi.nlm.nih.gov/16120773/)
    → **surface γ2 falls by about 47% one hour into SE.** `KENDO`/`KREC`/`KDEGR` were fitted to this single point, and the 30-minute value (0.69) is a result. CAL, PARAM

26. Goodkin HP, Yeh JL, Kapur J. **Status epilepticus increases the intracellular accumulation of GABA-A receptors.** J Neurosci 2005. — [PMID 15944379](https://pubmed.ncbi.nlm.nih.gov/15944379/)
    → the basis for the existence of the `RENDO` compartment (the endosomal pool)

27. **Goodkin HP, et al. Subunit-specific trafficking of GABA(A) receptors during status epilepticus.** J Neurosci 2008. — [PMID 18322097](https://pubmed.ncbi.nlm.nih.gov/18322097/)
    → **β2/3- and γ2-containing receptors internalise; δ-containing ones do not.** Why `RSYN` falls in the model and `REXTRA` does not. §3 vs §4

28. Terunuma M, et al. **Deficits in phosphorylation of GABA(A) receptors by intimately associated protein kinase C activity underlie compromised synaptic inhibition during status epilepticus.** J Neurosci 2008. — [PMID 18184780](https://pubmed.ncbi.nlm.nih.gov/18184780/)
    → the `T_ca` → `T_gephyrin` phosphorylation pathway

29. Kittler JT, et al. **Phospho-dependent binding of the clathrin AP2 adaptor complex to GABA-A receptors regulates the efficacy of inhibitory synaptic transmission.** PNAS 2005. — [PMID 16192353](https://pubmed.ncbi.nlm.nih.gov/16192353/)
    → the `T_ap2` node

30. Jacob TC, Moss SJ, Jurd R. **GABA(A) receptor trafficking and its role in the dynamic modulation of neuronal inhibition.** Nat Rev Neurosci 2008. — [PMID 18382465](https://pubmed.ncbi.nlm.nih.gov/18382465/)

31. Feng HJ, Mathews GC, Kao C, Macdonald RL. **Alterations of GABA-A receptor function and allosteric modulation during development of status epilepticus.** J Neurophysiol 2008. — [PMID 18216225](https://pubmed.ncbi.nlm.nih.gov/18216225/)
    → the basis for the **second** decreasing term, `FBZS` (the benzodiazepine-sensitive fraction)

32. Deeb TZ, Maguire J, Moss SJ. **Possible alterations in GABA-A receptor signaling that underlie benzodiazepine-resistant seizures.** Epilepsia 2012. — [PMID 23216581](https://pubmed.ncbi.nlm.nih.gov/23216581/)

33. Goodkin HP, Kapur J. **The impact of diazepam's discovery on the treatment and understanding of status epilepticus.** Epilepsia 2009. — [PMID 19674049](https://pubmed.ncbi.nlm.nih.gov/19674049/)

34. Joshi S, Kapur J. **Mechanisms of status epilepticus: α-amino-3-hydroxy-5-methyl-4-isoxazolepropionic acid receptor hypothesis.** Epilepsia 2018. — [PMID 30159880](https://pubmed.ncbi.nlm.nih.gov/30159880/)
    → the `AMPACP` (increase in Ca-permeable AMPA) state variable

---

## 4. **The pool that remains** — extrasynaptic δ-GABA-A and the neurosteroids

35. Farrant M, Nusser Z. **Variations on an inhibitory theme: phasic and tonic activation of GABA(A) receptors.** Nat Rev Neurosci 2005. — [PMID 15738957](https://pubmed.ncbi.nlm.nih.gov/15738957/)
    → the distinction between phasic (synaptic) and tonic (extrasynaptic) — the model's `WSYN`/`WEXT` decomposition

36. Rogawski MA, Loya CM, Reddy K, et al. **Neuroactive steroids for the treatment of status epilepticus.** Epilepsia 2013. — [PMID 24001085](https://pubmed.ncbi.nlm.nih.gov/24001085/)
    → δ receptors are benzodiazepine-insensitive and highly responsive to neurosteroids. `FSYNALLO = 0.40`

37. Reddy DS, Rogawski MA. **Neurosteroid replacement therapy for catamenial epilepsy.** Neurotherapeutics 2009. — [PMID 19332335](https://pubmed.ncbi.nlm.nih.gov/19332335/) → the `D_allowith` node

38. Zolkowska D, et al. **Intramuscular allopregnanolone and ganaxolone in a mouse model of treatment-resistant status epilepticus.** Epilepsia 2018. — [PMID 29453777](https://pubmed.ncbi.nlm.nih.gov/29453777/)

39. Rosenthal ES, et al. **Brexanolone as adjunctive therapy in super-refractory status epilepticus.** Ann Neurol 2017. — [PMID 28779545](https://pubmed.ncbi.nlm.nih.gov/28779545/)
    → the dose (86 µg/kg/h) and timing of scenario 18

40. Vaitkevicius H, et al. **Intravenous ganaxolone for the treatment of refractory status epilepticus: results from an open-label, dose-finding, phase 2 trial.** Epilepsia 2022. — [PMID 35748707](https://pubmed.ncbi.nlm.nih.gov/35748707/)

41. Singh RK, et al. **Intravenous ganaxolone in pediatric super-refractory status epilepticus.** Epilepsy Behav Rep 2022. — [PMID 36325100](https://pubmed.ncbi.nlm.nih.gov/36325100/)

42. Gasior M, et al. **Intravenous ganaxolone: pharmacokinetics, pharmacodynamics, safety, and tolerability in healthy adults.** Clin Pharmacol Drug Dev 2024. — [PMID 38231434](https://pubmed.ncbi.nlm.nih.gov/38231434/)

---

## 5. **Clock 2** — the arrival of synaptic NMDA receptors, and ketamine

43. **Naylor DE, et al. Rapid surface accumulation of NMDA receptors increases glutamatergic excitation during status epilepticus.** Neurobiol Dis 2013. — [PMID 23313318](https://pubmed.ncbi.nlm.nih.gov/23313318/)
    → **synaptic NR1 rises by about 38% one hour into SE.** `KEXO`/`KENDN`/`FEXO0` were fitted to this single point. CAL

44. Mazarati AM, Wasterlain CG. **N-methyl-D-aspartate receptor antagonists abolish the maintenance phase of self-sustaining status epilepticus in rat.** Neurosci Lett 1999. — [PMID 10327162](https://pubmed.ncbi.nlm.nih.gov/10327162/)
    → NMDA blockade does work in the maintenance phase = why the time sign of ketamine is positive

45. Rice AC, DeLorenzo RJ. **NMDA receptor activation during status epilepticus is required for the development of epilepsy.** Brain Res 1998. — [PMID 9519269](https://pubmed.ncbi.nlm.nih.gov/9519269/) → §19 `Z_epilep`

46. Fujikawa DG. **Neuroprotective effect of ketamine administered after status epilepticus onset.** Epilepsia 1995. — [PMID 7821277](https://pubmed.ncbi.nlm.nih.gov/7821277/)

47. Gaspard N, et al. **Intravenous ketamine for the treatment of refractory status epilepticus: a retrospective multicenter study.** Epilepsia 2013. — [PMID 23758557](https://pubmed.ncbi.nlm.nih.gov/23758557/)
    → the ketamine dose (a 1.5-3 mg/kg bolus, 1-10 mg/kg/h). PARAM

48. Rosati A, et al. **Efficacy and safety of ketamine in refractory status epilepticus in children.** Neurology 2012. — [PMID 23197747](https://pubmed.ncbi.nlm.nih.gov/23197747/)

49. Rosati A, De Masi S, Guerrini R. **Ketamine for refractory status epilepticus: a systematic review.** CNS Drugs 2018. — [PMID 30232735](https://pubmed.ncbi.nlm.nih.gov/30232735/)

50. Alkhachroum A, et al. **Ketamine to treat super-refractory status epilepticus.** Neurology 2020. — [PMID 32873691](https://pubmed.ncbi.nlm.nih.gov/32873691/)

51. Ilvento L, et al. **Ketamine in refractory convulsive status epilepticus in children avoids endotracheal intubation.** Epilepsy Behav 2015. — [PMID 26189786](https://pubmed.ncbi.nlm.nih.gov/26189786/)
    → why in the model ketamine alone **raises** MAP (`KETMAP`) and does not add to the respiratory burden

52. **Niquet J, et al. Simultaneous triple therapy for the treatment of status epilepticus.** Neurobiol Dis 2017. — [PMID 28461248](https://pubmed.ncbi.nlm.nih.gov/28461248/)
    → combination therapy aimed at **both** moving pools at once — the basis for scenarios 15 and 17

---

## 6. Chloride homeostasis — the route by which GABA becomes excitatory

53. Rivera C, et al. **BDNF-induced TrkB activation down-regulates the K⁺-Cl⁻ cotransporter KCC2 and impairs neuronal Cl⁻ extrusion.** J Cell Biol 2002. — [PMID 12473684](https://pubmed.ncbi.nlm.nih.gov/12473684/) → `C_bdnf` → `KCC2`

54. Pathak HR, et al. **Disrupted dentate granule cell chloride regulation enhances synaptic excitability during development of temporal lobe epilepsy.** J Neurosci 2007. — [PMID 18094240](https://pubmed.ncbi.nlm.nih.gov/18094240/)

55. Kaila K, et al. **Cation-chloride cotransporters in neuronal development, plasticity and disease.** Nat Rev Neurosci 2014. — [PMID 25234263](https://pubmed.ncbi.nlm.nih.gov/25234263/)

56. **Burman RJ, et al. Excitatory GABAergic signalling is associated with benzodiazepine resistance in status epilepticus.** Brain 2019. — [PMID 31553050](https://pubmed.ncbi.nlm.nih.gov/31553050/)
    → the fourth failure mode. The basis in the model for the feedback in which **the GABA agonist itself loads Cl⁻** (`FCLDRG`)

57. Pressler RM, et al. **Bumetanide for the treatment of seizures in newborn babies with hypoxic ischaemic encephalopathy (NEMO).** Lancet Neurol 2015. — [PMID 25765333](https://pubmed.ncbi.nlm.nih.gov/25765333/)
    → the basis for the clinical failure of NKCC1 targeting — stated at the `C_nkcc1` node

---

## 7. Endogenous terminators and promoters

58. Young D, Dragunow M. **Status epilepticus may be caused by loss of adenosine anticonvulsant mechanisms.** Neuroscience 1994. — [PMID 8152537](https://pubmed.ncbi.nlm.nih.gov/8152537/) → the `ADO` state variable

59. Aronica E, et al. **Glial adenosine kinase — a neuropathological marker of the epileptic brain.** Neurochem Int 2013. — [PMID 23385089](https://pubmed.ncbi.nlm.nih.gov/23385089/) → the rise in `ADK`

60. Mazarati AM, et al. **Galanin modulation of seizures and seizure modulation of hippocampal galanin in animal models of status epilepticus.** J Neurosci 1998. — [PMID 9822761](https://pubmed.ncbi.nlm.nih.gov/9822761/) → depletion of `NETPEP`

61. Mazarati AM, et al. **Modulation of hippocampal excitability and seizures by galanin.** J Neurosci 2000. — [PMID 10934278](https://pubmed.ncbi.nlm.nih.gov/10934278/)

62. Liu H, et al. **Substance P is expressed in hippocampal principal neurons during status epilepticus and plays a critical role in the maintenance of status epilepticus.** PNAS 1999. — [PMID 10220458](https://pubmed.ncbi.nlm.nih.gov/10220458/)
    → the inhibitory peptides are depleted while the excitatory ones rise = the monotonic fall of `NETPEP`

---

## 8. Neuroinflammation · the blood-brain barrier

63. Vezzani A, French J, Bartfai T, Baram TZ. **The role of inflammation in epilepsy.** Nat Rev Neurol 2011. — [PMID 21135885](https://pubmed.ncbi.nlm.nih.gov/21135885/)

64. Maroso M, et al. **Toll-like receptor 4 and high-mobility group box-1 are involved in ictogenesis and can be targeted to reduce seizures.** Nat Med 2010. — [PMID 20348922](https://pubmed.ncbi.nlm.nih.gov/20348922/) → `I_hmgb1`

65. **Balosso S, et al. A novel non-transcriptional pathway mediates the proconvulsive effects of interleukin-1β.** Brain 2008. — [PMID 18952671](https://pubmed.ncbi.nlm.nih.gov/18952671/)
    → IL-1β → Src → NR2B Tyr1472 phosphorylation → increased NMDA current. The model's `GIL` term

66. Friedman A, Kaufer D, Heinemann U. **Blood-brain barrier breakdown-inducing astrocytic transformation: novel targets for the prevention of epilepsy.** Epilepsy Res 2009. — [PMID 19362806](https://pubmed.ncbi.nlm.nih.gov/19362806/) → albumin-TGF-β → `Z_epilep`

67. Marchi N, et al. **Seizure-promoting effect of blood-brain barrier disruption.** Epilepsia 2007. — [PMID 17319915](https://pubmed.ncbi.nlm.nih.gov/17319915/)

68. Kenney-Jung DL, et al. **Febrile infection-related epilepsy syndrome treated with anakinra.** Ann Neurol 2016. — [PMID 27770579](https://pubmed.ncbi.nlm.nih.gov/27770579/) → scenario 20

69. Lai YC, et al. **Anakinra usage in febrile infection related epilepsy syndrome: an international cohort.** Ann Clin Transl Neurol 2020. — [PMID 33506622](https://pubmed.ncbi.nlm.nih.gov/33506622/)

70. Gaspard N, et al. **New-onset refractory status epilepticus: etiology, clinical features, and outcome.** Neurology 2015. — [PMID 26296517](https://pubmed.ncbi.nlm.nih.gov/26296517/)

71. Sculier C, Gaspard N. **New onset refractory status epilepticus (NORSE).** Seizure 2019. — [PMID 30482654](https://pubmed.ncbi.nlm.nih.gov/30482654/)

---

## 9. Excitotoxic neuronal injury — the second integrator

72. **Meldrum BS, Brierley JB. Prolonged epileptic seizures in primates: ischemic cell change and its relation to ictal physiological events.** Arch Neurol 1973. — [PMID 4629379](https://pubmed.ncbi.nlm.nih.gov/4629379/)
    → injury occurs even under ventilation and maintained blood pressure = why `INJURY` has to be a state variable **separate** from `SEIZ`

73. Meldrum BS, Horton RW. **Physiology of status epilepticus in primates.** Arch Neurol 1973. — [PMID 4629380](https://pubmed.ncbi.nlm.nih.gov/4629380/) → the §10 → §11 transition

74. Wasterlain CG, Fujikawa DG, Penix L, Sankar R. **Pathophysiological mechanisms of brain damage from status epilepticus.** Epilepsia 1993. — [PMID 8385002](https://pubmed.ncbi.nlm.nih.gov/8385002/)

75. Fujikawa DG. **The temporal evolution of neuronal damage from pilocarpine-induced status epilepticus.** Brain Res 1996. — [PMID 8828581](https://pubmed.ncbi.nlm.nih.gov/8828581/) → the `NEURLOSS` saturation curve

76. DeGiorgio CM, et al. **Serum neuron-specific enolase in the major subtypes of status epilepticus.** Neurology 1999. — [PMID 10078721](https://pubmed.ncbi.nlm.nih.gov/10078721/) → the `X_nse` node

---

## 10. Prognosis · outcome measures

77. Rossetti AO, Logroscino G, Bromfield EB. **Status Epilepticus Severity Score (STESS).** J Neurol 2008. — [PMID 18769858](https://pubmed.ncbi.nlm.nih.gov/18769858/)

78. Leitinger M, et al. **Epidemiology-based mortality score in status epilepticus (EMSE).** Neurocrit Care 2015. — [PMID 25412806](https://pubmed.ncbi.nlm.nih.gov/25412806/)

79. Leitinger M, et al. **Predicting outcome of status epilepticus.** Epilepsy Behav 2015. — [PMID 26071999](https://pubmed.ncbi.nlm.nih.gov/26071999/)

80. Giovannini G, et al. **Mortality, morbidity and refractoriness prediction in status epilepticus: comparison of STESS and EMSE scores.** Seizure 2017. — [PMID 28226274](https://pubmed.ncbi.nlm.nih.gov/28226274/)

81. Hocker S, Tatum WO, LaRoche S, Freeman WD. **Refractory and super-refractory status epilepticus — an update.** Curr Neurol Neurosci Rep 2014. — [PMID 24760477](https://pubmed.ncbi.nlm.nih.gov/24760477/)

82. Pujar SS, et al. **Long-term prognosis after childhood convulsive status epilepticus: a prospective cohort study.** Lancet Child Adolesc Health 2018. — [PMID 30169233](https://pubmed.ncbi.nlm.nih.gov/30169233/) → §19 `Z_epilep`

---

## 11. First-line treatment trials

83. **Treiman DM, et al. A comparison of four treatments for generalized convulsive status epilepticus (VA Cooperative Study).** N Engl J Med 1998. — [PMID 9738086](https://pubmed.ncbi.nlm.nih.gov/9738086/)
    → **lorazepam 64.9% in overt SE, 7.7-24.2% in subtle SE.** The same drug, a different window. Dissociation 1 of §18 of the model, and the `MOTG`/`MOTDRG` parameters. CAL

84. **Alldredge BK, et al. A comparison of lorazepam, diazepam, and placebo for the treatment of out-of-hospital status epilepticus (PHTSE).** N Engl J Med 2001. — [PMID 11547716](https://pubmed.ncbi.nlm.nih.gov/11547716/)
    → lorazepam 59.1% / diazepam 42.6% / placebo 21.1%; **respiratory complications were 22.5% on placebo against 10.6% on lorazepam** — the seizure depresses respiration more than the drug does. The basis for `KRESPS > KRESPD`. CAL

85. **Silbergleit R, et al. Intramuscular versus intravenous therapy for prehospital status epilepticus (RAMPART).** N Engl J Med 2012. — [PMID 22335736](https://pubmed.ncbi.nlm.nih.gov/22335736/)
    → **IM midazolam 73.4% vs IV lorazepam 63.4%.** The model has not a single parameter dedicated to RAMPART — it reproduces the result purely by the arithmetic that the absorption delay is smaller than the administration delay (model 80.0 vs 78.8, direction right, margin underestimated). CAL

86. Welch RD, et al. **Intramuscular midazolam versus intravenous lorazepam for the prehospital treatment of status epilepticus in the pediatric population.** Epilepsia 2015. — [PMID 25597369](https://pubmed.ncbi.nlm.nih.gov/25597369/)

87. Leppik IE, et al. **Double-blind study of lorazepam and diazepam in status epilepticus.** JAMA 1983. — [PMID 6131148](https://pubmed.ncbi.nlm.nih.gov/6131148/) → the redistribution and early recurrence of diazepam

88. McMullan J, et al. **Midazolam versus diazepam for the treatment of status epilepticus in children and young adults: a meta-analysis.** Acad Emerg Med 2010. — [PMID 20624136](https://pubmed.ncbi.nlm.nih.gov/20624136/)

89. **Sathe AG, et al. Underdosing of benzodiazepines in patients with status epilepticus enrolled in the Established Status Epilepticus Treatment Trial.** Acad Emerg Med 2019. — [PMID 31161706](https://pubmed.ncbi.nlm.nih.gov/31161706/)
    → underdosing is the rule, not the exception. The basis for scenarios 05 and 06

90. Sathe AG, et al. **Patterns of benzodiazepine underdosing in the Established Status Epilepticus Treatment Trial.** Epilepsia 2021. — [PMID 33567109](https://pubmed.ncbi.nlm.nih.gov/33567109/)

---

## 12. Second- and third-line treatment trials

91. **Kapur J, et al. Randomized trial of three anticonvulsant medications for status epilepticus (ESETT).** N Engl J Med 2019. — [PMID 31774955](https://pubmed.ncbi.nlm.nih.gov/31774955/)
    → **levetiracetam 47% / fosphenytoin 45% / valproate 46%.** The model's `SLEV`, `SPHT`, and `SVPA` (and their corresponding EC50s) are the only efficacy parameters fitted, and they were fitted to this **tie** alone (model 46.6 / 50.3 / 47.2). CAL

92. Chamberlain JM, et al. **Efficacy of levetiracetam, fosphenytoin, and valproate for established status epilepticus by age group (ESETT).** Lancet 2020. — [PMID 32203691](https://pubmed.ncbi.nlm.nih.gov/32203691/)

93. Navarro V, et al. **Prehospital treatment with levetiracetam plus clonazepam or placebo plus clonazepam in status epilepticus (SAMUKeppra).** Lancet Neurol 2016. — [PMID 26627366](https://pubmed.ncbi.nlm.nih.gov/26627366/)
    → adding pre-hospital levetiracetam was negative — consistent with the model's prediction that the benefit of a second-line agent is confined to those who failed the first line

94. Lyttle MD, et al. **Levetiracetam versus phenytoin for second-line treatment of paediatric convulsive status epilepticus (EcLiPSE).** Lancet 2019. — [PMID 31005385](https://pubmed.ncbi.nlm.nih.gov/31005385/)

95. Dalziel SR, et al. **Levetiracetam versus phenytoin for second-line treatment of convulsive status epilepticus in children (ConSEPT).** Lancet 2019. — [PMID 31005386](https://pubmed.ncbi.nlm.nih.gov/31005386/)

96. Claassen J, Hirsch LJ, Emerson RG, Mayer SA. **Treatment of refractory status epilepticus with pentobarbital, propofol, or midazolam: a systematic review.** Epilepsia 2002. — [PMID 11903460](https://pubmed.ncbi.nlm.nih.gov/11903460/) → §14

97. Rossetti AO, et al. **A randomized trial for the treatment of refractory status epilepticus.** Neurocrit Care 2011. — [PMID 20878265](https://pubmed.ncbi.nlm.nih.gov/20878265/)

98. Prabhakar H, et al. **Propofol versus thiopental sodium for the treatment of refractory status epilepticus.** Cochrane Database Syst Rev 2017. — [PMID 28155226](https://pubmed.ncbi.nlm.nih.gov/28155226/)

99. Legriel S, et al. **Hypothermia for neuroprotection in convulsive status epilepticus (HYBERNATUS).** N Engl J Med 2016. — [PMID 28002714](https://pubmed.ncbi.nlm.nih.gov/28002714/)
    → EEG seizures fell but functional outcome did not improve — why the `COOL` flag acts on `INJURY` only and not on `EDRIVE`

---

## 13. Pharmacokinetics (sources of the PK parameters)

100. Greenblatt DJ, et al. **Pharmacokinetic comparison of sublingual lorazepam with intravenous, intramuscular, and oral lorazepam.** J Pharm Sci 1982. — [PMID 6121043](https://pubmed.ncbi.nlm.nih.gov/6121043/) → `VLZP`, `CLLZP`

101. Ramsay RE, DeToledo J. **Intravenous administration of fosphenytoin: options for the management of seizures.** Neurology 1996. — [PMID 8649609](https://pubmed.ncbi.nlm.nih.gov/8649609/) → `KCONV`, the infusion rate limit

102. Clements JA, Nimmo WS. **Pharmacokinetics and analgesic effect of ketamine in man.** Br J Anaesth 1981. — [PMID 7459184](https://pubmed.ncbi.nlm.nih.gov/7459184/) → `V1KET`/`V2KET`/`CLKET`

103. Patsalos PN, et al. **Serum protein binding of 25 antiepileptic drugs in a routine clinical setting: a comparison of free non-protein-bound concentrations.** Epilepsia 2017. — [PMID 28542801](https://pubmed.ncbi.nlm.nih.gov/28542801/)
     → `FUPHT0 = 0.10`, and the saturable binding of valproate (`FUVPA0`/`FUVPAM`/`FUVPA50`)

104. Patsalos PN, et al. **Therapeutic drug monitoring of antiepileptic drugs in epilepsy: a 2018 update.** Ther Drug Monit 2018. — [PMID 29957667](https://pubmed.ncbi.nlm.nih.gov/29957667/)

---

## 14. Drug resistance — the transporter hypothesis vs the target hypothesis

105. Löscher W, Potschka H. **Drug resistance in brain diseases and the role of drug efflux transporters.** Nat Rev Neurosci 2005. — [PMID 16025095](https://pubmed.ncbi.nlm.nih.gov/16025095/)

106. Remy S, Beck H. **Molecular and cellular mechanisms of pharmacoresistance in epilepsy.** Brain 2006. — [PMID 16317026](https://pubmed.ncbi.nlm.nih.gov/16317026/)
     → the distinction between the two hypotheses. The model's `BPRATIO` output separates them

107. Bankstahl JP, Löscher W. **Resistance to antiepileptic drugs and expression of P-glycoprotein in two rat models of status epilepticus.** Epilepsy Res 2008. — [PMID 18760905](https://pubmed.ncbi.nlm.nih.gov/18760905/) → the `PGP` state variable

108. Bauer B, et al. **Seizure-induced up-regulation of P-glycoprotein at the blood-brain barrier through glutamate and cyclooxygenase-2 signaling.** Mol Pharmacol 2008. — [PMID 18094072](https://pubmed.ncbi.nlm.nih.gov/18094072/) → the `I_cox2` → `R_pgp` edge

109. Zibell G, et al. **Prevention of seizure-induced up-regulation of endothelial P-glycoprotein by COX-2 inhibition.** Neuropharmacology 2009. — [PMID 19371577](https://pubmed.ncbi.nlm.nih.gov/19371577/)

110. Potschka H, Fedrowitz M, Löscher W. **P-glycoprotein-mediated efflux of phenobarbital, lamotrigine, and felbamate at the blood-brain barrier.** Neurosci Lett 2002. — [PMID 12113905](https://pubmed.ncbi.nlm.nih.gov/12113905/)
     → the basis for the P-gp sensitivity of each drug (`FPGPBZD` ≪ `FPGPPHT`)

---

## 15. Non-convulsive SE and the EEG — the demonstration of dissociation 1

111. Treiman DM, Walton NY, Kendrick C. **A progressive sequence of electroencephalographic changes during generalized convulsive status epilepticus.** Epilepsy Res 1990. — [PMID 2303022](https://pubmed.ncbi.nlm.nih.gov/2303022/)
     → the EEG progresses even after the motor manifestations disappear — the direct basis for the decay of `MOTG`

112. DeLorenzo RJ, et al. **Persistent nonconvulsive status epilepticus after the control of convulsive status epilepticus.** Epilepsia 1998. — [PMID 9701373](https://pubmed.ncbi.nlm.nih.gov/9701373/)
     → EEG seizures persisted in 48% of the patients judged clinically to have stopped

113. Hirsch LJ, et al. **American Clinical Neurophysiology Society's Standardized Critical Care EEG Terminology: 2021 version.** J Clin Neurophysiol 2021. — [PMID 33475321](https://pubmed.ncbi.nlm.nih.gov/33475321/)

114. Claassen J, et al. **Detection of electrographic seizures with continuous EEG monitoring in critically ill patients.** Neurology 2004. — [PMID 15159471](https://pubmed.ncbi.nlm.nih.gov/15159471/)

115. Leitinger M, et al. **Salzburg consensus criteria for non-convulsive status epilepticus — approach to clinical application.** Epilepsy Behav 2015. — [PMID 26092326](https://pubmed.ncbi.nlm.nih.gov/26092326/)

116. Trinka E, Leitinger M. **Which EEG patterns in coma are nonconvulsive status epilepticus?** Epilepsy Behav 2015. — [PMID 26148985](https://pubmed.ncbi.nlm.nih.gov/26148985/)

117. Shneker BF, Fountain NB. **Assessment of acute morbidity and mortality in nonconvulsive status epilepticus.** Neurology 2003. — [PMID 14581666](https://pubmed.ncbi.nlm.nih.gov/14581666/)

---

## 16. The delay chain — the most powerful variable in the model, and the only purely operational one

118. **Gaínza-Lein M, et al. Association of time to treatment with short-term outcomes for pediatric patients with refractory convulsive status epilepticus.** JAMA Neurol 2018. — [PMID 29356811](https://pubmed.ncbi.nlm.nih.gov/29356811/)
     → the delay to the first dose is associated with death and with seizure duration independently — the clinical counterpart of `Q_cliff`

119. Sánchez Fernández I, et al. **Factors associated with treatment delays in pediatric refractory convulsive status epilepticus.** Neurology 2018. — [PMID 29643084](https://pubmed.ncbi.nlm.nih.gov/29643084/)

120. Kämppi L, Ritvanen J, Mustonen H, Soinila S. **Analysis of the delay components in the treatment of status epilepticus.** Neurocrit Care 2013. — [PMID 23817962](https://pubmed.ncbi.nlm.nih.gov/23817962/)
     → the measured components of the delay chain of §17

121. Kämppi L, et al. **The essence of the first 2.5 h in the treatment of generalized convulsive status epilepticus.** Seizure 2018. — [PMID 29306214](https://pubmed.ncbi.nlm.nih.gov/29306214/)

122. Guterman EL, et al. **Prehospital midazolam use and outcomes among patients with out-of-hospital status epilepticus.** Neurology 2020. — [PMID 32943481](https://pubmed.ncbi.nlm.nih.gov/32943481/)
     → real-world data showing that pre-hospital benzodiazepine is in fact given to only a minority

123. Vasquez A, et al. **Hospital emergency treatment of convulsive status epilepticus: comparison of pathways from ten pediatric research centers.** Pediatr Neurol 2018. — [PMID 30075875](https://pubmed.ncbi.nlm.nih.gov/30075875/)

---

## 17. Systemic physiology — the compensated and decompensated phases

124. Simon RP. **Physiologic consequences of status epilepticus.** Epilepsia 1985. — [PMID 3922751](https://pubmed.ncbi.nlm.nih.gov/3922751/)
     → `CATAMP`/`TAUCAT`, the failure of autoregulation, and the time course of hypoglycaemia and hyperthermia

125. Walton NY. **Systemic effects of generalized convulsive status epilepticus.** Epilepsia 1993. — [PMID 8462491](https://pubmed.ncbi.nlm.nih.gov/8462491/)

126. Manno EM, et al. **Cardiac pathology in status epilepticus.** Ann Neurol 2005. — [PMID 16240367](https://pubmed.ncbi.nlm.nih.gov/16240367/) → the `S2_takot` node

---

## Four points at which this model can be refuted (falsifiable predictions)

These are outputs of the model that **have not yet been verified** against the literature. If any one of them is wrong, the structure is wrong.

| # | Prediction | The value the model gives | How to test it |
|---|------|-----------------|-----------|
| 1 | A **vertical asymptote** in the benzodiazepine concentration required | between 45 and 50 minutes in the default patient | measure prospectively the point at which the first-line response rate collapses as a function of the time of administration |
| 2 | The cliff for the second-line agent **exists but is about 40 minutes later** | 25 minutes first line vs 62 minutes second line | response rates stratified by the time of second-line administration in an ESETT-type trial |
| 3 | The **crossing time** of the recoverable benzodiazepine pool and the removable NMDA share | 280 minutes | a trial randomising the time at which ketamine is added |
| 4 | The **difference in cumulative excitotoxic burden** between termination at 20 minutes and at 60 minutes | 9.3% vs 42.7% neuronal loss | follow-up of NSE/NfL and hippocampal volume by time of termination |

---

## Disclaimer

This model is a **qualitative to semi-quantitative QSP model for educational and research purposes**.
It was built on the published literature but has not been independently verified or certified, and
**it must not be used in actual clinical decision-making, in prescribing, or in regulatory submission.**
