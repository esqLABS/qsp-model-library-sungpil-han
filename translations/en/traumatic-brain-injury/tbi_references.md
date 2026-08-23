# Severe Traumatic Brain Injury — References
## References for the sTBI QSP model

The literature is arranged to follow the structure of the model. Each section states which equation or
parameter of the model came from which reference. **132 papers** in total.

> Citation convention: PubMed links are given by PMID; entries without a PMID are given by DOI.

---

## 1. Intracranial pressure-volume relationship and compliance (Monro-Kellie)

*Model correspondence: `PVI`, `Cic = PVI/(ICP·ln10)`, the `V_csf_min` buffering gate, the R2 reserve-capacity result*

1. Monro A. *Observations on the Structure and Functions of the Nervous System.* Edinburgh: Creech & Johnson, 1783. — the original source.
2. Kellie G. An account of the appearances observed in the dissection of two of three individuals... *Trans Med Chir Soc Edinb.* 1824;1:84-169.
3. Marmarou A, Shulman K, LaMorgese J. Compartmental analysis of compliance and outflow resistance of the cerebrospinal fluid system. *J Neurosurg.* 1975;43(5):523-34. [PMID 1181384](https://pubmed.ncbi.nlm.nih.gov/1181384/) — **the original source of PVI. Where PVI ≈ 26 mL comes from.**
4. Marmarou A, Shulman K, Rosende RM. A nonlinear analysis of the cerebrospinal fluid system and intracranial pressure dynamics. *J Neurosurg.* 1978;48(3):332-44. [PMID 632857](https://pubmed.ncbi.nlm.nih.gov/632857/) — the exponential P–V curve.
5. Löfgren J, von Essen C, Zwetnow NN. The pressure-volume curve of the cerebrospinal fluid space in dogs. *Acta Neurol Scand.* 1973;49(5):557-74. [PMID 4763787](https://pubmed.ncbi.nlm.nih.gov/4763787/)
6. Avezaat CJ, van Eijndhoven JH, Wyper DJ. Cerebrospinal fluid pulse pressure and intracranial volume-pressure relationships. *J Neurol Neurosurg Psychiatry.* 1979;42(8):687-700. [PMID 490174](https://pubmed.ncbi.nlm.nih.gov/490174/)
7. Czosnyka M, Pickard JD. Monitoring and interpretation of intracranial pressure. *J Neurol Neurosurg Psychiatry.* 2004;75(6):813-21. [PMID 15145991](https://pubmed.ncbi.nlm.nih.gov/15145991/) — comprehensive review.
8. Czosnyka M, Guazzo E, Whitehouse M, et al. Significance of intracranial pressure waveform analysis after head injury. *Acta Neurochir.* 1996;138(5):531-41. [PMID 8800328](https://pubmed.ncbi.nlm.nih.gov/8800328/) — the RAP index.
9. Kim DJ, Czosnyka Z, Keong N, et al. Index of cerebrospinal compensatory reserve in hydrocephalus. *Neurosurgery.* 2009;64(3):494-501. [PMID 19240611](https://pubmed.ncbi.nlm.nih.gov/19240611/)
10. Wilson MH. Monro-Kellie 2.0: The dynamic vascular and venous pathophysiological components of intracranial pressure. *J Cereb Blood Flow Metab.* 2016;36(8):1338-50. [PMID 27174995](https://pubmed.ncbi.nlm.nih.gov/27174995/) — **basis for the venous buffering compartment (`Vv_max`).**
11. Mokri B. The Monro-Kellie hypothesis: applications in CSF volume depletion. *Neurology.* 2001;56(12):1746-8. [PMID 11425944](https://pubmed.ncbi.nlm.nih.gov/11425944/)
12. Piper IR, Miller JD, Dearden NM, et al. Systems analysis of cerebrovascular pressure transmission. *J Neurosurg.* 1990;73(6):871-80. [PMID 2230970](https://pubmed.ncbi.nlm.nih.gov/2230970/)

---

## 2. Cerebral autoregulation and the pressure reactivity index (PRx, CPPopt)

*Model correspondence: `x_aut`, the compliance sigmoid `Ca(x)`, `G_aut`, the `PRx` bisection solution, R1·R4·R5*

13. Lassen NA. Cerebral blood flow and oxygen consumption in man. *Physiol Rev.* 1959;39(2):183-238. [PMID 13645234](https://pubmed.ncbi.nlm.nih.gov/13645234/) — **the original source of the Lassen curve.**
14. Ursino M, Lodi CA. A simple mathematical model of the interaction between intracranial pressure and cerebral hemodynamics. *J Appl Physiol.* 1997;82(4):1256-69. [PMID 9104864](https://pubmed.ncbi.nlm.nih.gov/9104864/) — **the structural foundation of the haemodynamic core of this model.**
15. Ursino M, Lodi CA. Interaction among autoregulation, CO2 reactivity, and intracranial pressure: a mathematical model. *Am J Physiol.* 1998;274(5):H1715-28. [PMID 9612384](https://pubmed.ncbi.nlm.nih.gov/9612384/) — **the interaction between CO2 and autoregulation.**
16. Ursino M. A mathematical study of human intracranial hydrodynamics. Part 1: The cerebrospinal fluid pulse pressure. *Ann Biomed Eng.* 1988;16(4):379-401. [PMID 3177984](https://pubmed.ncbi.nlm.nih.gov/3177984/)
17. Czosnyka M, Smielewski P, Kirkpatrick P, et al. Continuous assessment of the cerebral vasomotor reactivity in head injury. *Neurosurgery.* 1997;41(1):11-9. [PMID 9218290](https://pubmed.ncbi.nlm.nih.gov/9218290/) — **the original source of PRx.**
18. Steiner LA, Czosnyka M, Piechnik SK, et al. Continuous monitoring of cerebrovascular pressure reactivity allows determination of optimal cerebral perfusion pressure. *Crit Care Med.* 2002;30(4):733-8. [PMID 11940737](https://pubmed.ncbi.nlm.nih.gov/11940737/) — **CPPopt.**
19. Aries MJ, Czosnyka M, Budohoski KP, et al. Continuous determination of optimal cerebral perfusion pressure in traumatic brain injury. *Crit Care Med.* 2012;40(8):2456-63. [PMID 22622398](https://pubmed.ncbi.nlm.nih.gov/22622398/)
20. Needham E, McFadyen C, Newcombe V, et al. Cerebral perfusion pressure targets individualized to pressure-reactivity index in moderate to severe TBI: a systematic review. *J Neurotrauma.* 2017;34(5):963-70. [PMID 27246184](https://pubmed.ncbi.nlm.nih.gov/27246184/)
21. Beqiri E, Smielewski P, Robba C, et al. Feasibility of individualised severe traumatic brain injury management using an automated assessment of optimal cerebral perfusion pressure: COGiTATE phase II study. *BMJ Open.* 2019;9(9):e030727. [PMID 31509426](https://pubmed.ncbi.nlm.nih.gov/31509426/)
22. Tas J, Beqiri E, van Kaam RC, et al. Targeting autoregulation-guided cerebral perfusion pressure after traumatic brain injury (COGiTATE): a feasibility randomized controlled clinical trial. *J Neurotrauma.* 2021;38(20):2790-800. [PMID 34407385](https://pubmed.ncbi.nlm.nih.gov/34407385/)
23. Drummond JC. The lower limit of autoregulation: time to revise our thinking? *Anesthesiology.* 1997;86(6):1431-3. [PMID 9197320](https://pubmed.ncbi.nlm.nih.gov/9197320/)
24. Rosner MJ, Rosner SD, Johnson AH. Cerebral perfusion pressure: management protocol and clinical results. *J Neurosurg.* 1995;83(6):949-62. [PMID 7490638](https://pubmed.ncbi.nlm.nih.gov/7490638/) — **Rosner's CPP-augmentation theory (R4).**
25. Rosner MJ, Becker DP. Origin and evolution of plateau waves. Experimental observations and a theoretical model. *J Neurosurg.* 1984;60(2):312-24. [PMID 6693959](https://pubmed.ncbi.nlm.nih.gov/6693959/) — **the vasodilatory cascade and plateau waves (R3).**
26. Lundberg N. Continuous recording and control of ventricular fluid pressure in neurosurgical practice. *Acta Psychiatr Scand Suppl.* 1960;36(149):1-193. [PMID 13764297](https://pubmed.ncbi.nlm.nih.gov/13764297/) — **the original source of the A wave (plateau wave).**
27. Castellani G, Zweifel C, Kim DJ, et al. Plateau waves in head injured patients requiring neurocritical care. *Neurocrit Care.* 2009;11(2):143-50. [PMID 19479205](https://pubmed.ncbi.nlm.nih.gov/19479205/)
28. Dias C, Maia I, Cerejo A, et al. Pressures, flow, and brain oxygenation during plateau waves of intracranial pressure. *Neurocrit Care.* 2014;21(1):124-32. [PMID 24343564](https://pubmed.ncbi.nlm.nih.gov/24343564/)
29. Nordström CH, Reinstrup P, Xu W, et al. Assessment of the lower limit for cerebral perfusion pressure in severe head injuries by bedside monitoring of regional energy metabolism. *Anesthesiology.* 2003;98(4):809-14. [PMID 12657839](https://pubmed.ncbi.nlm.nih.gov/12657839/)
30. Grände PO. The "Lund Concept" for the treatment of severe head trauma. *Intensive Care Med.* 2006;32(10):1475-84. [PMID 16896859](https://pubmed.ncbi.nlm.nih.gov/16896859/) — **the Lund concept (the R19 comparator arm).**
31. Zeiler FA, Ercole A, Cabeleira M, et al. Comparison of performance of different optimal cerebral perfusion pressure parameters for outcome prediction in adult TBI: CENTER-TBI. *J Neurotrauma.* 2019;36(10):1505-17. [PMID 30384809](https://pubmed.ncbi.nlm.nih.gov/30384809/)
32. Rivera-Lara L, Zorrilla-Vaca A, Geocadin RG, et al. Cerebral autoregulation-oriented therapy at the bedside: a comprehensive review. *Anesthesiology.* 2017;126(6):1187-99. [PMID 28361741](https://pubmed.ncbi.nlm.nih.gov/28361741/)

---

## 3. CO₂ reactivity and its bicarbonate adaptation

*Model correspondence: `f_co2`, `PaCO2_eff = PaCO2·24/HCO3`, `tau_hco3`, `k_co2_cal`, R8*

33. Kety SS, Schmidt CF. The effects of altered arterial tensions of carbon dioxide and oxygen on cerebral blood flow and cerebral oxygen consumption of normal young men. *J Clin Invest.* 1948;27(4):484-92. [PMID 16695569](https://pubmed.ncbi.nlm.nih.gov/16695569/) — **the original source of CO2 reactivity.**
34. Kontos HA, Raper AJ, Patterson JL. Analysis of vasoactivity of local pH, PCO2 and bicarbonate on pial vessels. *Stroke.* 1977;8(3):358-60. [PMID 17497](https://pubmed.ncbi.nlm.nih.gov/17497/) — **vessels respond to perivascular pH, not to PaCO2 — the basis for `PaCO2_eff` in this model.**
35. Muizelaar JP, Marmarou A, Ward JD, et al. Adverse effects of prolonged hyperventilation in patients with severe head injury: a randomized clinical trial. *J Neurosurg.* 1991;75(5):731-9. [PMID 1919695](https://pubmed.ncbi.nlm.nih.gov/1919695/) — **the harm of sustained hyperventilation. The clinical counterpart of R8.**
36. Muizelaar JP, van der Poel HG, Li ZC, et al. Pial arteriolar vessel diameter and CO2 reactivity during prolonged hyperventilation in the rabbit. *J Neurosurg.* 1988;69(6):923-7. [PMID 3142974](https://pubmed.ncbi.nlm.nih.gov/3142974/) — **the mechanism of the loss of effect (bicarbonate adaptation).**
37. Raichle ME, Posner JB, Plum F. Cerebral blood flow during and after hyperventilation. *Arch Neurol.* 1970;23(5):394-403. [PMID 5471647](https://pubmed.ncbi.nlm.nih.gov/5471647/)
38. Coles JP, Fryer TD, Coleman MR, et al. Hyperventilation following head injury: effect on ischemic burden and cerebral oxidative metabolism. *Crit Care Med.* 2007;35(2):568-78. [PMID 17205016](https://pubmed.ncbi.nlm.nih.gov/17205016/) — **the oxygen debt of hyperventilation (the "interest rate" of R8).**
39. Coles JP, Minhas PS, Fryer TD, et al. Effect of hyperventilation on cerebral blood flow in traumatic head injury. *Crit Care Med.* 2002;30(9):1950-9. [PMID 12352026](https://pubmed.ncbi.nlm.nih.gov/12352026/)
40. Marion DW, Puccio A, Wisniewski SR, et al. Effect of hyperventilation on extracellular concentrations of glutamate, lactate, pyruvate, and local cerebral blood flow in patients with severe traumatic brain injury. *Crit Care Med.* 2002;30(12):2619-25. [PMID 12483051](https://pubmed.ncbi.nlm.nih.gov/12483051/)
41. Godoy DA, Seifi A, Chi G, et al. The intracranial compartmental syndrome: a proposed model for acute brain injury monitoring and management. *Crit Care.* 2023;27(1):137. [PMID 37038236](https://pubmed.ncbi.nlm.nih.gov/37038236/)
42. Brandi G, Stocchetti N, Pagnamenta A, et al. Cerebral metabolism is not affected by moderate hyperventilation in patients with traumatic brain injury. *Crit Care.* 2019;23(1):45. [PMID 30760282](https://pubmed.ncbi.nlm.nih.gov/30760282/)

---

## 4. Cerebral metabolism, oxygen delivery, brain tissue oxygen tension (DO₂, PbtO₂, SjvO₂)

*Model correspondence: `CMRO2n`, `Q10`, the `softmin` supply limit, `PbtO2`, `SjvO2`, R14*

43. Rosenthal G, Hemphill JC, Sorani M, et al. Brain tissue oxygen tension is more indicative of oxygen diffusion than oxygen delivery and metabolism in patients with traumatic brain injury. *Crit Care Med.* 2008;36(6):1917-24. [PMID 18496363](https://pubmed.ncbi.nlm.nih.gov/18496363/) — **basis for the `Dfac` diffusion factor.**
44. Maloney-Wilensky E, Gracias V, Itkin A, et al. Brain tissue oxygen and outcome after severe traumatic brain injury: a systematic review. *Crit Care Med.* 2009;37(6):2057-63. [PMID 19384213](https://pubmed.ncbi.nlm.nih.gov/19384213/)
45. Okonkwo DO, Shutter LA, Moore C, et al. Brain oxygen optimization in severe traumatic brain injury phase-II (BOOST-II): a phase II randomized trial. *Crit Care Med.* 2017;45(11):1907-14. [PMID 28841632](https://pubmed.ncbi.nlm.nih.gov/28841632/) — **BOOST-2.**
46. Bernard F, Barsan W, Diaz-Arrastia R, et al. Brain Oxygen Optimization in Severe Traumatic Brain Injury (BOOST-3): a multicentre, randomised, blinded-endpoint, comparative effectiveness study. *BMJ Open.* 2022;12(3):e060188. [PMID 35296490](https://pubmed.ncbi.nlm.nih.gov/35296490/) — **BOOST-3 (the basis for R14).**
47. Payen JF, Launey Y, Chabanne R, et al. Intracranial pressure monitoring with and without brain tissue oxygen pressure monitoring for severe traumatic brain injury in France (OXY-TC): an open-label, randomised controlled superiority trial. *Lancet Neurol.* 2023;22(11):1005-14. [PMID 37863590](https://pubmed.ncbi.nlm.nih.gov/37863590/)
48. Robertson CS, Narayan RK, Gokaslan ZL, et al. Cerebral arteriovenous oxygen difference as an estimate of cerebral blood flow in comatose patients. *J Neurosurg.* 1989;70(2):222-30. [PMID 2913221](https://pubmed.ncbi.nlm.nih.gov/2913221/)
49. Gopinath SP, Robertson CS, Contant CF, et al. Jugular venous desaturation and outcome after head injury. *J Neurol Neurosurg Psychiatry.* 1994;57(6):717-23. [PMID 8006653](https://pubmed.ncbi.nlm.nih.gov/8006653/) — **the prognostic significance of SjvO2.**
50. Coles JP, Fryer TD, Smielewski P, et al. Incidence and mechanisms of cerebral ischemia in early clinical head injury. *J Cereb Blood Flow Metab.* 2004;24(2):202-11. [PMID 14747747](https://pubmed.ncbi.nlm.nih.gov/14747747/) — **PET evidence that global indices miss regional ischaemia.**
51. Menon DK, Coles JP, Gupta AK, et al. Diffusion limited oxygen delivery following head injury. *Crit Care Med.* 2004;32(6):1384-90. [PMID 15187523](https://pubmed.ncbi.nlm.nih.gov/15187523/) — **microvascular diffusion limitation (`Dfac`, `pen_res`).**
52. Severinghaus JW. Simple, accurate equations for human blood O2 dissociation computations. *J Appl Physiol.* 1979;46(3):599-602. [PMID 35496](https://pubmed.ncbi.nlm.nih.gov/35496/) — **source of the inverse ODC equation.**
53. Vespa P, Bergsneider M, Hattori N, et al. Metabolic crisis without brain ischemia is common after traumatic brain injury: a combined microdialysis and positron emission tomography study. *J Cereb Blood Flow Metab.* 2005;25(6):763-74. [PMID 15716852](https://pubmed.ncbi.nlm.nih.gov/15716852/) — **metabolic crisis without ischaemia = a mitochondrial-type rise in LPR.**
54. Timofeev I, Carpenter KL, Nortje J, et al. Cerebral extracellular chemistry and outcome following traumatic brain injury: a microdialysis study of 223 patients. *Brain.* 2011;134(2):484-94. [PMID 21247930](https://pubmed.ncbi.nlm.nih.gov/21247930/) — **basis for the LPR threshold of 40.**
55. Hutchinson PJ, Jalloh I, Helmy A, et al. Consensus statement from the 2014 International Microdialysis Forum. *Intensive Care Med.* 2015;41(9):1517-28. [PMID 26194024](https://pubmed.ncbi.nlm.nih.gov/26194024/)
56. Vespa PM, McArthur D, O'Phelan K, et al. Persistently low extracellular glucose correlates with poor outcome 6 months after human traumatic brain injury. *J Cereb Blood Flow Metab.* 2003;23(7):865-77. [PMID 12843790](https://pubmed.ncbi.nlm.nih.gov/12843790/)
57. Oddo M, Schmidt JM, Carrera E, et al. Impact of tight glycemic control on cerebral glucose metabolism after severe brain injury: a microdialysis study. *Crit Care Med.* 2008;36(12):3233-8. [PMID 18936695](https://pubmed.ncbi.nlm.nih.gov/18936695/)
58. Michenfelder JD, Milde JH. The relationship among canine brain temperature, metabolism, and function during hypothermia. *Anesthesiology.* 1991;75(1):130-6. [PMID 2064037](https://pubmed.ncbi.nlm.nih.gov/2064037/) — **where Q10 ≈ 2.3 comes from.**

---

## 5. Cerebral oedema: cytotoxic and vasogenic (reflection coefficient, AQP4, SUR1)

*Model correspondence: the two-compartment Starling equation, `sig_prot_open`, `k_idio`, `Lh_inj`, R6*

59. Klatzo I. Presidental address. Neuropathological aspects of brain edema. *J Neuropathol Exp Neurol.* 1967;26(1):1-14. [PMID 5336776](https://pubmed.ncbi.nlm.nih.gov/5336776/) — **the original source of the cytotoxic/vasogenic distinction.**
60. Fenstermacher JD, Johnson JA. Filtration and reflection coefficients of the rabbit blood-brain barrier. *Am J Physiol.* 1966;211(2):341-6. [PMID 5921084](https://pubmed.ncbi.nlm.nih.gov/5921084/) — **measurement of the reflection coefficient σ ≈ 1. A core parameter of this model.**
61. Zlokovic BV. The blood-brain barrier in health and chronic neurodegenerative disorders. *Neuron.* 2008;57(2):178-201. [PMID 18215617](https://pubmed.ncbi.nlm.nih.gov/18215617/)
62. Stokum JA, Gerzanich V, Simard JM. Molecular pathophysiology of cerebral edema. *J Cereb Blood Flow Metab.* 2016;36(3):513-38. [PMID 26661240](https://pubmed.ncbi.nlm.nih.gov/26661240/) — comprehensive review.
63. Simard JM, Kent TA, Chen M, et al. Brain oedema in focal ischaemia: molecular pathophysiology and theoretical implications. *Lancet Neurol.* 2007;6(3):258-68. [PMID 17303532](https://pubmed.ncbi.nlm.nih.gov/17303532/)
64. Simard JM, Woo SK, Schwartzbauer GT, Gerzanich V. Sulfonylurea receptor 1 in central nervous system injury: a focused review. *J Cereb Blood Flow Metab.* 2012;32(9):1699-717. [PMID 22714048](https://pubmed.ncbi.nlm.nih.gov/22714048/) — **SUR1-TRPM4, the glibenclamide target.**
65. Jha RM, Bell J, Citerio G, et al. Role of sulfonylurea receptor 1 and glibenclamide in traumatic brain injury: a review of the evidence. *Int J Mol Sci.* 2020;21(2):409. [PMID 31936452](https://pubmed.ncbi.nlm.nih.gov/31936452/)
66. Manley GT, Fujimura M, Ma T, et al. Aquaporin-4 deletion in mice reduces brain edema after acute water intoxication and ischemic stroke. *Nat Med.* 2000;6(2):159-63. [PMID 10655103](https://pubmed.ncbi.nlm.nih.gov/10655103/) — **AQP4.**
67. Papadopoulos MC, Verkman AS. Aquaporin water channels in the nervous system. *Nat Rev Neurosci.* 2013;14(4):265-77. [PMID 23481483](https://pubmed.ncbi.nlm.nih.gov/23481483/)
68. Iliff JJ, Wang M, Liao Y, et al. A paravascular pathway facilitates CSF flow through the brain parenchyma and the clearance of interstitial solutes, including amyloid β. *Sci Transl Med.* 2012;4(147):147ra111. [PMID 22896675](https://pubmed.ncbi.nlm.nih.gov/22896675/) — **the glymphatic system (`k_glym`).**
69. Iliff JJ, Chen MJ, Plog BA, et al. Impairment of glymphatic pathway function promotes tau pathology after traumatic brain injury. *J Neurosci.* 2014;34(49):16180-93. [PMID 25471560](https://pubmed.ncbi.nlm.nih.gov/25471560/) — **reduced glymphatic clearance in the injured region (`k_glym_inj` < `k_glym_int`).**
70. Marmarou A, Signoretti S, Fatouros PP, et al. Predominance of cellular edema in traumatic brain swelling in patients with severe head injuries. *J Neurosurg.* 2006;104(5):720-30. [PMID 16703876](https://pubmed.ncbi.nlm.nih.gov/16703876/)
71. Unterberg AW, Stover J, Kress B, Kiening KL. Edema and brain trauma. *Neuroscience.* 2004;129(4):1021-9. [PMID 15561417](https://pubmed.ncbi.nlm.nih.gov/15561417/)
72. Donkin JJ, Vink R. Mechanisms of cerebral edema in traumatic brain injury: therapeutic developments. *Curr Opin Neurol.* 2010;23(3):293-9. [PMID 20168229](https://pubmed.ncbi.nlm.nih.gov/20168229/)

---

## 6. Osmotherapy: mannitol, hypertonic saline, rebound

*Model correspondence: `Mann_c/p/inj`, `Na_ecf`, `Osm_p`, the `k_mann_bbb` rebound mechanism, R6·R7·R9*

73. Wise BL, Chater N. The value of hypertonic mannitol solution in decreasing brain mass and lowering cerebro-spinal-fluid pressure. *J Neurosurg.* 1962;19:1038-43. [PMID 13998001](https://pubmed.ncbi.nlm.nih.gov/13998001/)
74. Kaufmann AM, Cardoso ER. Aggravation of vasogenic cerebral edema by multiple-dose mannitol. *J Neurosurg.* 1992;77(4):584-9. [PMID 1527618](https://pubmed.ncbi.nlm.nih.gov/1527618/) — **the experimental basis for mannitol rebound. The core of R7.**
75. Node Y, Nakazawa S. Clinical study of mannitol and glycerol on raised intracranial pressure and on their rebound phenomenon. *Adv Neurol.* 1990;52:359-63. [PMID 2109910](https://pubmed.ncbi.nlm.nih.gov/2109910/)
76. Diringer MN, Zazulia AR. Osmotic therapy: fact and fiction. *Neurocrit Care.* 2004;1(2):219-33. [PMID 16174920](https://pubmed.ncbi.nlm.nih.gov/16174920/)
77. Diringer MN, Scalfani MT, Zazulia AR, et al. Effect of mannitol on cerebral blood volume in patients with head injury. *Neurosurgery.* 2012;70(5):1215-9. [PMID 22089756](https://pubmed.ncbi.nlm.nih.gov/22089756/)
78. Cottenceau V, Masson F, Mahamid E, et al. Comparison of effects of equiosmolar doses of mannitol and hypertonic saline on cerebral blood flow and metabolism in traumatic brain injury. *J Neurotrauma.* 2011;28(10):2003-12. [PMID 21787184](https://pubmed.ncbi.nlm.nih.gov/21787184/) — **an equiosmolar comparison. The clinical counterpart of R6.**
79. Mangat HS, Wu X, Gerber LM, et al. Hypertonic saline is superior to mannitol for the combined effect on intracranial pressure and cerebral perfusion pressure burdens in patients with severe traumatic brain injury. *Neurosurgery.* 2020;86(2):221-30. [PMID 30924498](https://pubmed.ncbi.nlm.nih.gov/30924498/)
80. Roquilly A, Moyer JD, Huet O, et al. Effect of continuous infusion of hypertonic saline vs standard care on 6-month neurological outcomes in patients with traumatic brain injury: the COBI randomized clinical trial. *JAMA.* 2021;325(20):2056-66. [PMID 34032829](https://pubmed.ncbi.nlm.nih.gov/34032829/) — **continuous hypertonic saline did not improve outcome.**
81. Cooper DJ, Myles PS, McDermott FT, et al. Prehospital hypertonic saline resuscitation of patients with hypotension and severe traumatic brain injury: a randomized controlled trial. *JAMA.* 2004;291(11):1350-7. [PMID 15026402](https://pubmed.ncbi.nlm.nih.gov/15026402/)
82. Bullock MR. Mannitol and other diuretics in severe neurotrauma. *New Horiz.* 1995;3(3):448-52. [PMID 7496753](https://pubmed.ncbi.nlm.nih.gov/7496753/)
83. Gondim FA, Aiyagari V, Shackleford A, Diringer MN. Osmolality not predictive of mannitol-induced acute renal insufficiency. *J Neurosurg.* 2005;103(3):444-7. [PMID 16235675](https://pubmed.ncbi.nlm.nih.gov/16235675/)
84. Visweswaran P, Massin EK, Dubose TD. Mannitol-induced acute renal failure. *J Am Soc Nephrol.* 1997;8(6):1028-33. [PMID 9189872](https://pubmed.ncbi.nlm.nih.gov/9189872/) — **nephrotoxicity at an osmolar gap > 55 (`Osm_aki`).**
85. García-Morales EJ, Cariappa R, Parvin CA, et al. Osmole gap in neurologic-neurosurgical intensive care unit: its normal value, calculation, and relationship with mannitol serum concentrations. *Crit Care Med.* 2004;32(4):986-91. [PMID 15071390](https://pubmed.ncbi.nlm.nih.gov/15071390/)
86. Rickard AC, Smith JE, Newell P, et al. Salt or sugar for your injured brain? A meta-analysis of randomised controlled trials of mannitol versus hypertonic sodium solutions. *Emerg Med J.* 2014;31(8):679-83. [PMID 23811861](https://pubmed.ncbi.nlm.nih.gov/23811861/)
87. Chen H, Song Z, Dennis JA. Hypertonic saline versus other intracranial pressure-lowering agents for people with acute traumatic brain injury. *Cochrane Database Syst Rev.* 2020;1:CD010904. [PMID 31978260](https://pubmed.ncbi.nlm.nih.gov/31978260/)
88. Ropper AH. Hyperosmolar therapy for raised intracranial pressure. *N Engl J Med.* 2012;367(8):746-52. [PMID 22913684](https://pubmed.ncbi.nlm.nih.gov/22913684/)

---

## 7. Excitotoxicity, mitochondria, spreading depolarisation

*Model correspondence: `Glu`, `Ca_i`, `MitoD`, `ROS`, `SD_rate`, the cell-death term `kill`*

89. Faden AI, Demediuk P, Panter SS, Vink R. The role of excitatory amino acids and NMDA receptors in traumatic brain injury. *Science.* 1989;244(4906):798-800. [PMID 2567056](https://pubmed.ncbi.nlm.nih.gov/2567056/)
90. Bullock R, Zauner A, Woodward JJ, et al. Factors affecting excitatory amino acid release following severe human head injury. *J Neurosurg.* 1998;89(4):507-18. [PMID 9761042](https://pubmed.ncbi.nlm.nih.gov/9761042/) — **measured glutamate values (2 → 50-100 µM).**
91. Nilsson P, Hillered L, Pontén U, Ungerstedt U. Changes in cortical extracellular levels of energy-related metabolites and amino acids following concussive brain injury in rats. *J Cereb Blood Flow Metab.* 1990;10(5):631-7. [PMID 2384536](https://pubmed.ncbi.nlm.nih.gov/2384536/)
92. Sullivan PG, Rabchevsky AG, Waldmeier PC, Springer JE. Mitochondrial permeability transition in CNS trauma: cause or effect of neuronal cell death? *J Neurosci Res.* 2005;79(1-2):231-9. [PMID 15573402](https://pubmed.ncbi.nlm.nih.gov/15573402/)
93. Verweij BH, Muizelaar JP, Vinas FC, et al. Impaired cerebral mitochondrial function after traumatic brain injury in humans. *J Neurosurg.* 2000;93(5):815-20. [PMID 11059663](https://pubmed.ncbi.nlm.nih.gov/11059663/)
94. Astrup J, Siesjö BK, Symon L. Thresholds in cerebral ischemia — the ischemic penumbra. *Stroke.* 1981;12(6):723-5. [PMID 6272455](https://pubmed.ncbi.nlm.nih.gov/6272455/) — **two thresholds, electrical failure and ion-pump failure — the dual structure of the `kill` term in this model.**
95. Hossmann KA. Viability thresholds and the penumbra of focal ischemia. *Ann Neurol.* 1994;36(4):557-65. [PMID 7944288](https://pubmed.ncbi.nlm.nih.gov/7944288/)
96. Jones TH, Morawetz RB, Crowell RM, et al. Thresholds of focal cerebral ischemia in awake monkeys. *J Neurosurg.* 1981;54(6):773-82. [PMID 7241187](https://pubmed.ncbi.nlm.nih.gov/7241187/) — **CBF thresholds (18, 10-12 mL/100g/min).**
97. Hartings JA, Bullock MR, Okonkwo DO, et al. Spreading depolarisations and outcome after traumatic brain injury: a prospective observational study. *Lancet Neurol.* 2011;10(12):1058-64. [PMID 22014379](https://pubmed.ncbi.nlm.nih.gov/22014379/) — **COSBID.**
98. Dreier JP. The role of spreading depression, spreading depolarization and spreading ischemia in neurological disease. *Nat Med.* 2011;17(4):439-47. [PMID 21475241](https://pubmed.ncbi.nlm.nih.gov/21475241/) — **the inverse haemodynamic response (`SDINV`).**
99. Hinzman JM, Wilson JA, Mazzeo AT, et al. Excitotoxicity and metabolic crisis are associated with spreading depolarizations in severe traumatic brain injury patients. *J Neurotrauma.* 2016;33(19):1775-83. [PMID 26586606](https://pubmed.ncbi.nlm.nih.gov/26586606/)
100. Hartings JA, Shuttleworth CW, Kirov SA, et al. The continuum of spreading depolarizations in acute cortical lesion development. *J Cereb Blood Flow Metab.* 2017;37(5):1571-94. [PMID 27328690](https://pubmed.ncbi.nlm.nih.gov/27328690/)

---

## 8. Neuroinflammation and the blood-brain barrier (biphasic BBB opening)

*Model correspondence: `Micro`, `Cyto`, `Neut`, `MMP9`, `BBB_mech` + `BBB_infl`*

101. Baskaya MK, Rao AM, Doğan A, et al. The biphasic opening of the blood-brain barrier in the cortex and hippocampus after traumatic brain injury in rats. *Neurosci Lett.* 1997;226(1):33-6. [PMID 9153635](https://pubmed.ncbi.nlm.nih.gov/9153635/) — **biphasic barrier opening — the two state variables of this model.**
102. Shlosberg D, Benifla M, Kaufer D, Friedman A. Blood-brain barrier breakdown as a therapeutic target in traumatic brain injury. *Nat Rev Neurol.* 2010;6(7):393-403. [PMID 20551947](https://pubmed.ncbi.nlm.nih.gov/20551947/)
103. Corps KN, Roth TL, McGavern DB. Inflammation and neuroprotection in traumatic brain injury. *JAMA Neurol.* 2015;72(3):355-62. [PMID 25599342](https://pubmed.ncbi.nlm.nih.gov/25599342/)
104. Helmy A, Carpenter KL, Menon DK, et al. The cytokine response to human traumatic brain injury: temporal profiles and evidence for cerebral parenchymal production. *J Cereb Blood Flow Metab.* 2011;31(2):658-70. [PMID 20717122](https://pubmed.ncbi.nlm.nih.gov/20717122/)
105. Simon DW, McGeachy MJ, Bayır H, et al. The far-reaching scope of neuroinflammation after traumatic brain injury. *Nat Rev Neurol.* 2017;13(3):171-91. [PMID 28186177](https://pubmed.ncbi.nlm.nih.gov/28186177/)
106. Vajtr D, Benada O, Kukacka J, et al. Correlation of ultrastructural changes of endothelial cells and astrocytes occurring during blood brain barrier damage after traumatic brain injury with biochemical markers. *Physiol Res.* 2009;58(2):263-8. [PMID 18380546](https://pubmed.ncbi.nlm.nih.gov/18380546/)
107. Wang X, Jung J, Asahi M, et al. Effects of matrix metalloproteinase-9 gene knock-out on morphological and motor outcomes after traumatic brain injury. *J Neurosci.* 2000;20(18):7037-42. [PMID 10995849](https://pubmed.ncbi.nlm.nih.gov/10995849/)
108. Thelin EP, Tajsic T, Zeiler FA, et al. Monitoring the neuroinflammatory response following acute brain injury. *Front Neurol.* 2017;8:351. [PMID 28775710](https://pubmed.ncbi.nlm.nih.gov/28775710/)

---

## 9. Haematoma expansion, coagulopathy, tranexamic acid (TXA)

*Model correspondence: `V_hem`, `Fib`, `TXA_c/p`, `txa_IC50`, R15*

109. Maegele M, Schöchl H, Menovsky T, et al. Coagulopathy and haemorrhagic progression in traumatic brain injury: advances in mechanisms, diagnosis, and management. *Lancet Neurol.* 2017;16(8):630-47. [PMID 28721927](https://pubmed.ncbi.nlm.nih.gov/28721927/)
110. CRASH-3 trial collaborators. Effects of tranexamic acid on death, disability, vascular occlusive events and other morbidities in patients with acute traumatic brain injury (CRASH-3): a randomised, placebo-controlled trial. *Lancet.* 2019;394(10210):1713-23. [PMID 31623894](https://pubmed.ncbi.nlm.nih.gov/31623894/) — **the clinical counterpart of R15.**
111. Gayet-Ageron A, Prieto-Merino D, Ker K, et al. Effect of treatment delay on the effectiveness and safety of antifibrinolytics in acute severe haemorrhage: a meta-analysis of individual patient-level data from 40 138 bleeding patients. *Lancet.* 2018;391(10116):125-32. [PMID 29126600](https://pubmed.ncbi.nlm.nih.gov/29126600/) — **the curve of declining effect with treatment delay.**
112. Narayan RK, Maas AI, Servadei F, et al. Progression of traumatic intracerebral hemorrhage: a prospective observational study. *J Neurotrauma.* 2008;25(6):629-39. [PMID 18491950](https://pubmed.ncbi.nlm.nih.gov/18491950/)
113. Kurland D, Hong C, Aarabi B, et al. Hemorrhagic progression of a contusion after traumatic brain injury: a review. *J Neurotrauma.* 2012;29(1):19-31. [PMID 21988198](https://pubmed.ncbi.nlm.nih.gov/21988198/)
114. Zhang J, Zhang F, Dong JF. Coagulopathy induced by traumatic brain injury: systemic manifestation of a localized injury. *Blood.* 2018;131(18):2001-6. [PMID 29507078](https://pubmed.ncbi.nlm.nih.gov/29507078/) — **release of brain tissue factor.**

---

## 10. Guidelines, interventional trials, tiered therapy

*Model correspondence: the `Tiered` closed-loop protocol, R10·R11·R12*

115. Carney N, Totten AM, O'Reilly C, et al. Guidelines for the Management of Severe Traumatic Brain Injury, Fourth Edition. *Neurosurgery.* 2017;80(1):6-15. [PMID 27654000](https://pubmed.ncbi.nlm.nih.gov/27654000/) — **BTF 4th edition. Source of ICP 22 mmHg and CPP 60-70.**
116. Hawryluk GWJ, Aguilera S, Buki A, et al. A management algorithm for patients with intracranial pressure monitoring: the Seattle International Severe Traumatic Brain Injury Consensus Conference (SIBICC). *Intensive Care Med.* 2019;45(12):1783-94. [PMID 31659383](https://pubmed.ncbi.nlm.nih.gov/31659383/) — **source of the R10 tiered algorithm.**
117. Chesnut RM, Temkin N, Carney N, et al. A trial of intracranial-pressure monitoring in traumatic brain injury (BEST:TRIP). *N Engl J Med.* 2012;367(26):2471-81. [PMID 23234472](https://pubmed.ncbi.nlm.nih.gov/23234472/) — **monitoring is not in itself a treatment.**
118. Cooper DJ, Rosenfeld JV, Murray L, et al. Decompressive craniectomy in diffuse traumatic brain injury (DECRA). *N Engl J Med.* 2011;364(16):1493-502. [PMID 21434843](https://pubmed.ncbi.nlm.nih.gov/21434843/) — **early surgery worsened outcome (R11).**
119. Hutchinson PJ, Kolias AG, Timofeev IS, et al. Trial of decompressive craniectomy for traumatic intracranial hypertension (RESCUEicp). *N Engl J Med.* 2016;375(12):1119-30. [PMID 27602507](https://pubmed.ncbi.nlm.nih.gov/27602507/) — **last-resort surgery: fewer deaths, more severe disability (R11).**
120. Andrews PJ, Sinclair HL, Rodriguez A, et al. Hypothermia for intracranial hypertension after traumatic brain injury (Eurotherm3235). *N Engl J Med.* 2015;373(25):2403-12. [PMID 26444221](https://pubmed.ncbi.nlm.nih.gov/26444221/) — **ICP came down and the outcome got worse (the admitted limitation of R12).**
121. Cooper DJ, Nichol AD, Bailey M, et al. Effect of early sustained prophylactic hypothermia on neurologic outcomes among patients with severe traumatic brain injury: the POLAR randomized clinical trial. *JAMA.* 2018;320(21):2211-20. [PMID 30357266](https://pubmed.ncbi.nlm.nih.gov/30357266/)
122. Clifton GL, Valadka A, Zygun D, et al. Very early hypothermia induction in patients with severe brain injury (NABIS: H II). *Lancet Neurol.* 2011;10(2):131-9. [PMID 21169065](https://pubmed.ncbi.nlm.nih.gov/21169065/)
123. Roberts I, Sydenham E. Barbiturates for acute traumatic brain injury. *Cochrane Database Syst Rev.* 2012;12:CD000033. [PMID 23235573](https://pubmed.ncbi.nlm.nih.gov/23235573/)
124. Robertson CS, Valadka AB, Hannay HJ, et al. Prevention of secondary ischemic insults after severe head injury. *Crit Care Med.* 1999;27(10):2086-95. [PMID 10548187](https://pubmed.ncbi.nlm.nih.gov/10548187/) — **a CPP 70 target increased ARDS.**
125. Chesnut RM, Marshall LF, Klauber MR, et al. The role of secondary brain injury in determining outcome from severe head injury. *J Trauma.* 1993;34(2):216-22. [PMID 8459458](https://pubmed.ncbi.nlm.nih.gov/8459458/) — **a single hypotensive episode doubles mortality (R13).**
126. Spaite DW, Hu C, Bobrow BJ, et al. Mortality and prehospital blood pressure in patients with major traumatic brain injury: implications for the hypotension threshold. *JAMA Surg.* 2017;152(4):360-8. [PMID 27926759](https://pubmed.ncbi.nlm.nih.gov/27926759/) — **shows that the relationship is continuous rather than a threshold.**
127. Güiza F, Depreitere B, Piper I, et al. Visualizing the pressure and time burden of intracranial hypertension in adult and paediatric traumatic brain injury. *Intensive Care Med.* 2015;41(6):1067-76. [PMID 25894624](https://pubmed.ncbi.nlm.nih.gov/25894624/) — **source of the pressure-time burden (dose) concept.**
128. Vik A, Nag T, Fredriksli OA, et al. Relationship of "dose" of intracranial hypertension to outcome in severe traumatic brain injury. *J Neurosurg.* 2008;109(4):678-84. [PMID 18826355](https://pubmed.ncbi.nlm.nih.gov/18826355/)

---

## 11. Biomarkers and prognostic models (IMPACT, CRASH, GOS-E)

*Model correspondence: the `GFAP`/`UCHL1`/`NfL`/`S100B` half-lives, `gose_unfavourable()`*

129. Bazarian JJ, Biberthaler P, Welch RD, et al. Serum GFAP and UCH-L1 for prediction of absence of intracranial injuries on head CT (ALERT-TBI). *Lancet Neurol.* 2018;17(9):782-9. [PMID 30054151](https://pubmed.ncbi.nlm.nih.gov/30054151/) — **GFAP/UCH-L1 half-lives and time profiles.**
130. Papa L, Brophy GM, Welch RD, et al. Time course and diagnostic accuracy of glial and neuronal blood biomarkers GFAP and UCH-L1 in a large cohort of trauma patients. *JAMA Neurol.* 2016;73(5):551-60. [PMID 27018834](https://pubmed.ncbi.nlm.nih.gov/27018834/) — **UCH-L1 peaks at 8 h, GFAP at 20-24 h — the comparator for R16.**
131. Shahim P, Politis A, van der Merwe A, et al. Neurofilament light as a biomarker in traumatic brain injury. *Neurology.* 2020;95(6):e610-22. [PMID 32641529](https://pubmed.ncbi.nlm.nih.gov/32641529/)
132. Steyerberg EW, Mushkudiani N, Perel P, et al. Predicting outcome after traumatic brain injury: development and international validation of prognostic scores based on admission characteristics (IMPACT). *PLoS Med.* 2008;5(8):e165. [PMID 18684008](https://pubmed.ncbi.nlm.nih.gov/18684008/) — **the structural prototype of the prognostic logistic in this model.**

---

## Appendix A — Source cross-reference table for the model parameters

| Parameter | Value | Source |
|---|---|---|
| `PVI` pressure-volume index | 26 mL | Marmarou 1975 [3] |
| `Can` arteriolar compliance | 0.20 mL/mmHg | Ursino & Lodi 1997 [14] |
| `G_aut` autoregulation gain | 3.0 | Ursino & Lodi 1997 [14], recalibrated in this model to the Lassen plateau |
| `tau_aut` | 0.333 min (20 s) | Ursino & Lodi 1997 [14] |
| CO2 reactivity | 3.2 %/mmHg (30-40) | Kety & Schmidt 1948 [33] |
| `tau_hco3` bicarbonate adaptation | 360 min | Muizelaar 1988 [36] |
| `Q10` | 2.3 | Michenfelder & Milde 1991 [58] |
| `CMRO2n` | 3.3 mL O₂/100g/min | Kety & Schmidt 1948 [33] |
| Reflection coefficient σ (intact BBB) | ≈ 0.97 | Fenstermacher & Johnson 1966 [60] |
| CBF ischaemia threshold | 18 / 10-12 mL/100g/min | Astrup 1981 [94], Jones 1981 [96] |
| LPR abnormality threshold | 40 | Timofeev 2011 [54] |
| PbtO₂ hypoxia threshold | 20 mmHg | Maloney-Wilensky 2009 [44] |
| ICP treatment threshold | 22 mmHg | Carney 2017 (BTF 4th edition) [115] |
| CPP target | 60-70 mmHg | Carney 2017 [115] |
| Osmolar gap nephrotoxicity threshold | 55 mOsm/kg | Visweswaran 1997 [84] |
| GFAP / UCH-L1 half-life | ~24 h / ~7 h | Papa 2016 [130] |
| Tiered therapy algorithm | tier 0-3 | Hawryluk 2019 (SIBICC) [116] |
| Prognostic model structure | age·pupils·CT | Steyerberg 2008 (IMPACT) [132] |

## Appendix B — Literature the model fails to reproduce

For the sake of honesty, the points at which the model is **wrong** are set out explicitly.

| Reference | What was observed | What the model predicts | Why they differ |
|---|---|---|---|
| Eurotherm3235 [120], POLAR [121] | Hypothermia lowers ICP but worsens outcome | It both lowers ICP and improves outcome | The model contains no pneumonia, no immunosuppression, no coagulopathy, no shivering, and no arrhythmia. A model that holds only the mechanisms by which a treatment is believed to work has no choice but to predict that the treatment works. |
| DECRA [118] | Early decompression worsened outcome | It reproduces the small benefit of early surgery, but not the net harm | There are no surgical complications, no syndrome of the trephined, and no reoperation. |
| COBI [80] | Continuous hypertonic saline did not improve outcome | It predicts a reduction in ICP burden | The link this model itself makes — that a reduced ICP burden is an improved outcome — is what is under test. |

This table is a list of the model's defects and, at the same time, tells us what this model is a tool for.
A mechanistic model is a tool for computing **what happens if the mechanisms are right**,
not a tool for deciding **whether that treatment benefits the patient**.
