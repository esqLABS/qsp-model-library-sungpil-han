# Dementia with Lewy Bodies (DLB) — Reference List

Supporting literature list (for `dlb_qsp_model.dot`,
`dlb_mrgsolve_model.R`, `dlb_shiny_app_en.R`).

**Link convention.** Every entry below resolves through a **title-resolving PubMed
search URL**(`pubmed.ncbi.nlm.nih.gov/?term=...`). Fixed PMIDs are deliberately not
asserted, because in an automatically generated session a wrong PMID mixed in makes
tracing harder rather than easier. Each link lands directly on the paper concerned (or on the
result set in which it appears first).

*Every entry below links to a **title-resolving PubMed query**. Fixed PMIDs are
deliberately not asserted: a single mistyped accession number is worse than a
query that always resolves. Each link lands on the paper (or a result set with
it first).*

The `[model location]` marker after each entry indicates which part of the model that
reference supports — a parameter, an equation, or a calibration anchor.

---

## 1. Diagnostic criteria and clinical phenotype

1. McKeith IG, Boeve BF, Dickson DW, et al. **Diagnosis and management of dementia with Lewy bodies: Fourth consensus report of the DLB Consortium.** *Neurology* 2017. — The core clinical features (cognitive fluctuation · visual hallucinations · REM sleep behaviour disorder · parkinsonism) and the indicative biomarkers (DaTSCAN · MIBG · PSG). The basis for the whole of the model's clinical endpoint definitions. [`$INIT` clinical states, `$TABLE`]
   <https://pubmed.ncbi.nlm.nih.gov/?term=Diagnosis+and+management+of+dementia+with+Lewy+bodies+Fourth+consensus+report+of+the+DLB+Consortium>
2. McKeith IG, Ferman TJ, Thomas AJ, et al. **Research criteria for the diagnosis of prodromal dementia with Lewy bodies.** *Neurology* 2020. — Definition of the prodromal stage (MCI-LB · psychiatric-onset · delirium-onset). The conceptual basis for the "prodromal dosing" of scenario 19.
   <https://pubmed.ncbi.nlm.nih.gov/?term=Research+criteria+for+the+diagnosis+of+prodromal+dementia+with+Lewy+bodies>
3. Walker Z, Possin KL, Boeve BF, Aarsland D. **Lewy body dementias.** *Lancet* 2015. — The view of DLB and Parkinson's disease dementia (PDD) as one Lewy body disease spectrum. The basis for the design of the model's `PHENO` switch.
   <https://pubmed.ncbi.nlm.nih.gov/?term=Lewy+body+dementias+Walker+Possin+Boeve+Aarsland+Lancet>
4. Ferman TJ, Smith GE, Boeve BF, et al. **DLB fluctuations: specific features that reliably differentiate DLB from AD and normal aging.** *Neurology* 2004. — The four-item scale (Mayo Fluctuations Composite) that distinguishes cognitive fluctuation from AD. The basis on which the model treats `CAF` as "variance rather than mean". [`FLUCTGT`]
   <https://pubmed.ncbi.nlm.nih.gov/?term=DLB+fluctuations+specific+features+that+reliably+differentiate+DLB+from+AD+and+normal+aging>
5. Walker MP, Ayre GA, Cummings JL, et al. **The Clinician Assessment of Fluctuation and the One Day Fluctuation Assessment Scale.** *Br J Psychiatry* 2000. — The original paper for the CAF scale (0–16). The basis for the model's `CAF = 4 × FLUC` scaling.
   <https://pubmed.ncbi.nlm.nih.gov/?term=The+Clinician+Assessment+of+Fluctuation+and+the+One+Day+Fluctuation+Assessment+Scale>
6. Vann Jones SA, O'Brien JT. **The prevalence and incidence of dementia with Lewy bodies: a systematic review of population and clinical studies.** *Psychol Med* 2014. — Prevalence and incidence; male predominance.
   <https://pubmed.ncbi.nlm.nih.gov/?term=The+prevalence+and+incidence+of+dementia+with+Lewy+bodies+a+systematic+review+of+population+and+clinical+studies>
7. Kane JPM, Surendranathan A, Bentley A, et al. **Clinical prevalence of Lewy body dementia.** *Alzheimers Res Ther* 2018.
   <https://pubmed.ncbi.nlm.nih.gov/?term=Clinical+prevalence+of+Lewy+body+dementia+Kane+Surendranathan>
8. Emre M, Aarsland D, Brown R, et al. **Clinical diagnostic criteria for dementia associated with Parkinson's disease.** *Mov Disord* 2007. — The PDD criteria and the "one-year rule". The definition of the `PHENO = 1` arm.
   <https://pubmed.ncbi.nlm.nih.gov/?term=Clinical+diagnostic+criteria+for+dementia+associated+with+Parkinson+disease+Emre+Aarsland>

## 2. Neuropathology and staging

9. Braak H, Del Tredici K, Rüb U, et al. **Staging of brain pathology related to sporadic Parkinson's disease.** *Neurobiol Aging* 2003. — The ascending brainstem→limbic→neocortex staging. The model's three-compartment regional structure (`FIBB/FIBL/FIBN`) and its sequential gating (`UPTL`, `UPTN`).
   <https://pubmed.ncbi.nlm.nih.gov/?term=Staging+of+brain+pathology+related+to+sporadic+Parkinson+disease+Braak+Del+Tredici>
10. Attems J, Toledo JB, Walker L, et al. **Neuropathological consensus criteria for the evaluation of Lewy pathology in post-mortem brains.** *Acta Neuropathol* 2021. — Grading of brainstem-predominant / limbic / neocortical distribution.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Neuropathological+consensus+criteria+for+the+evaluation+of+Lewy+pathology+in+post-mortem+brains>
11. Beach TG, Adler CH, Lue L, et al. **Unified staging system for Lewy body disorders.** *Acta Neuropathol* 2009.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Unified+staging+system+for+Lewy+body+disorders+Beach+Adler+Lue>
12. Irwin DJ, Grossman M, Weintraub D, et al. **Neuropathological and genetic correlates of survival and dementia onset in synucleinopathies.** *Lancet Neurol* 2017. — Concomitant AD pathology shortens survival. The model's `B_CROSS`/`BNS` and `PTAU` → survival path.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Neuropathological+and+genetic+correlates+of+survival+and+dementia+onset+in+synucleinopathies>
13. Borghammer P, Van Den Berge N. **Brain-first versus gut-first Parkinson's disease: a hypothesis.** *J Parkinsons Dis* 2019. — The two entry routes `S_GUT` / `S_OLF`.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Brain-first+versus+gut-first+Parkinson+disease+a+hypothesis+Borghammer>
14. Halliday GM, Holton JL, Revesz T, Dickson DW. **Neuropathology underlying clinical variability in patients with synucleinopathies.** *Acta Neuropathol* 2011.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Neuropathology+underlying+clinical+variability+in+patients+with+synucleinopathies>

## 3. α-Synuclein biophysics · aggregation kinetics · strains

15. Cohen SIA, Linse S, Luheshi LM, et al. **Proliferation of amyloid-β42 aggregates occurs through a secondary nucleation mechanism.** *PNAS* 2013. — The quantitative framework for secondary nucleation (catalysis on the fibril surface). The form of the model's `KNUC2 · ASYNM · FIB` term itself.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Proliferation+of+amyloid-beta42+aggregates+occurs+through+a+secondary+nucleation+mechanism>
16. Buell AK, Galvagnion C, Gaspar R, et al. **Solution conditions determine the relative importance of nucleation and growth processes in α-synuclein aggregation.** *PNAS* 2014. — The relative contributions of primary and secondary nucleation for α-synuclein.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Solution+conditions+determine+the+relative+importance+of+nucleation+and+growth+processes+in+alpha-synuclein+aggregation>
17. Meisl G, Kirkegaard JB, Arosio P, et al. **Molecular mechanisms of protein aggregation from global fitting of kinetic models.** *Nat Protoc* 2016. — The standard form of aggregation kinetic models (including capacity-limiting terms). The model's `(1 − OLIG)` and `(1 − FIB)` saturation terms.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Molecular+mechanisms+of+protein+aggregation+from+global+fitting+of+kinetic+models>
18. Winner B, Jappelli R, Maji SK, et al. **In vivo demonstration that α-synuclein oligomers are toxic.** *PNAS* 2011. — Oligomers (not fibrils) are the toxic species. The basis on which the model couples neuronal loss to `OLIG` and propagation to `FIB`.
    <https://pubmed.ncbi.nlm.nih.gov/?term=In+vivo+demonstration+that+alpha-synuclein+oligomers+are+toxic>
19. Schweighauser M, Shi Y, Tarutani A, et al. **Structures of α-synuclein filaments from multiple system atrophy.** *Nature* 2020. — Disease-specific differences in filament structure (cryo-EM). The model's `STRAIN` / `A_STRAIN` parameters.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Structures+of+alpha-synuclein+filaments+from+multiple+system+atrophy>
20. Yang Y, Shi Y, Schweighauser M, et al. **Structures of α-synuclein filaments from human brains with Lewy pathology.** *Nature* 2022. — The fold specific to Lewy pathology.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Structures+of+alpha-synuclein+filaments+from+human+brains+with+Lewy+pathology>
21. Peng C, Gathagan RJ, Covell DJ, et al. **Cellular milieu imparts distinct pathological α-synuclein strains in α-synucleinopathies.** *Nature* 2018.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Cellular+milieu+imparts+distinct+pathological+alpha-synuclein+strains+in+alpha-synucleinopathies>
22. Fujiwara H, Hasegawa M, Dohmae N, et al. **α-Synuclein is phosphorylated in synucleinopathy lesions.** *Nat Cell Biol* 2002. — pS129 accounts for >90% of Lewy body α-synuclein. `A_PS129`.
    <https://pubmed.ncbi.nlm.nih.gov/?term=alpha-Synuclein+is+phosphorylated+in+synucleinopathy+lesions>
23. Shahmoradian SH, Lewis AJ, Genoud C, et al. **Lewy pathology in Parkinson's disease consists of crowded organelles and lipid membranes.** *Nat Neurosci* 2019. — Lewy bodies are assemblies of membranes and organelles rather than simple protein aggregates. The reason `A_LB` is treated as "an inert sink or an active reservoir".
    <https://pubmed.ncbi.nlm.nih.gov/?term=Lewy+pathology+in+Parkinson+disease+consists+of+crowded+organelles+and+lipid+membranes>
24. Volpicelli-Daley LA, Luk KC, Patel TP, et al. **Exogenous α-synuclein fibrils induce Lewy body pathology leading to synaptic dysfunction and neuron death.** *Neuron* 2011.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Exogenous+alpha-synuclein+fibrils+induce+Lewy+body+pathology+leading+to+synaptic+dysfunction+and+neuron+death>

## 4. Cell-to-cell propagation

25. Desplats P, Lee HJ, Bae EJ, et al. **Inclusion formation and neuronal cell death through neuron-to-neuron transmission of α-synuclein.** *PNAS* 2009. — The `S_EXO` → `S_ISF` → `S_UPT` path.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Inclusion+formation+and+neuronal+cell+death+through+neuron-to-neuron+transmission+of+alpha-synuclein>
26. Mao X, Ou MT, Karuppagounder SS, et al. **Pathological α-synuclein transmission initiated by binding lymphocyte-activation gene 3.** *Science* 2016. — LAG3 receptor-mediated uptake.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Pathological+alpha-synuclein+transmission+initiated+by+binding+lymphocyte-activation+gene+3>
27. Holmes BB, DeVos SL, Kfoury N, et al. **Heparan sulfate proteoglycans mediate internalization and propagation of specific proteopathic seeds.** *PNAS* 2013.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Heparan+sulfate+proteoglycans+mediate+internalization+and+propagation+of+specific+proteopathic+seeds>
28. Kim S, Kwon SH, Kam TI, et al. **Transneuronal propagation of pathologic α-synuclein from the gut to the brain models Parkinson's disease.** *Neuron* 2019. — The gastric and vagal route.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Transneuronal+propagation+of+pathologic+alpha-synuclein+from+the+gut+to+the+brain+models+Parkinson+disease>
29. Xie L, Kang H, Xu Q, et al. **Sleep drives metabolite clearance from the adult brain.** *Science* 2013. — Sleep-dependent glymphatic clearance. The model's `GLYMF = f(SLD)` — the feedback by which damage to the RBD circuitry worsens seed clearance.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Sleep+drives+metabolite+clearance+from+the+adult+brain>
30. Iliff JJ, Wang M, Liao Y, et al. **A paravascular pathway facilitates CSF flow through the brain parenchyma and the clearance of interstitial solutes.** *Sci Transl Med* 2012.
    <https://pubmed.ncbi.nlm.nih.gov/?term=A+paravascular+pathway+facilitates+CSF+flow+through+the+brain+parenchyma+and+the+clearance+of+interstitial+solutes>

## 5. GBA1 · lysosome · ambroxol

31. Mazzulli JR, Xu YH, Sun Y, et al. **Gaucher disease glucocerebrosidase and α-synuclein form a bidirectional pathogenic loop in synucleinopathies.** *Cell* 2011. — **The model's GBA1 bistable switch itself.** GCase↓ → GlcCer↑ → oligomer stabilisation → block of GCase trafficking → GCase↓. [`L_SWITCH`, `GCTGT`, `KOLDE`]
    <https://pubmed.ncbi.nlm.nih.gov/?term=Gaucher+disease+glucocerebrosidase+and+alpha-synuclein+form+a+bidirectional+pathogenic+loop+in+synucleinopathies>
32. Nalls MA, Duran R, Lopez G, et al. **A multicenter study of glucocerebrosidase mutations in dementia with Lewy bodies.** *JAMA Neurol* 2013. — The odds ratio for GBA1 variants in DLB ≈ 8, stronger than in Parkinson's disease. The `GBAF` parameter.
    <https://pubmed.ncbi.nlm.nih.gov/?term=A+multicenter+study+of+glucocerebrosidase+mutations+in+dementia+with+Lewy+bodies>
33. Sidransky E, Nalls MA, Aasly JO, et al. **Multicenter analysis of glucocerebrosidase mutations in Parkinson's disease.** *N Engl J Med* 2009.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Multicenter+analysis+of+glucocerebrosidase+mutations+in+Parkinson+disease+Sidransky+Nalls>
34. Migdalska-Richards A, Daly L, Bezard E, Schapira AHV. **Ambroxol effects in glucocerebrosidase and α-synuclein transgenic mice.** *Ann Neurol* 2016. — Ambroxol raises CNS GCase activity and lowers α-synuclein.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Ambroxol+effects+in+glucocerebrosidase+and+alpha-synuclein+transgenic+mice>
35. Mullin S, Smith L, Lee K, et al. **Ambroxol for the treatment of patients with Parkinson disease with and without glucocerebrosidase gene mutations: a non-randomized, noncontrolled trial.** *JAMA Neurol* 2020. — 1.26 g/day, leucocyte GCase activity about +35%, CSF penetration confirmed. The calibration anchor for the model's `EMAXAMB`/`EC50AMB`.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Ambroxol+for+the+treatment+of+patients+with+Parkinson+disease+with+and+without+glucocerebrosidase+gene+mutations>
36. Silveira CRA, MacKinley J, Coleman K, et al. **Ambroxol as a novel disease-modifying treatment for Parkinson's disease dementia: protocol for a single-centre, randomized, double-blind, placebo-controlled trial.** *BMC Neurol* 2019. — The design of an ambroxol trial targeting PDD.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Ambroxol+as+a+novel+disease-modifying+treatment+for+Parkinson+disease+dementia+protocol>
37. Cuervo AM, Stefanis L, Fredenburg R, et al. **Impaired degradation of mutant α-synuclein by chaperone-mediated autophagy.** *Science* 2004. — CMA/LAMP2A. `L_CMA`.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Impaired+degradation+of+mutant+alpha-synuclein+by+chaperone-mediated+autophagy>
38. Jinn S, Drolet RE, Cramer PE, et al. **TMEM175 deficiency impairs lysosomal and mitochondrial function and increases α-synuclein aggregation.** *PNAS* 2017. — `G_TMEM` → `L_ACID`.
    <https://pubmed.ncbi.nlm.nih.gov/?term=TMEM175+deficiency+impairs+lysosomal+and+mitochondrial+function+and+increases+alpha-synuclein+aggregation>
39. Zunke F, Moise AC, Belur NR, et al. **Reversible conformational conversion of α-synuclein into toxic assemblies by glucosylceramide.** *Neuron* 2018. — GlcCer stabilises the oligomeric conformation. `KGLC50`.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Reversible+conformational+conversion+of+alpha-synuclein+into+toxic+assemblies+by+glucosylceramide>
40. Peterschmitt MJ, Saiki H, Hatano T, et al. **Safety, pharmacokinetics, and pharmacodynamics of oral venglustat in Parkinson disease patients with a GBA mutation.** *Mov Disord* 2022. — A GCS inhibitor (upstream blockade). `Y_VENG`.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Safety+pharmacokinetics+and+pharmacodynamics+of+oral+venglustat+in+Parkinson+disease+patients+with+a+GBA+mutation>

## 6. Genetics beyond GBA1

41. Guerreiro R, Ross OA, Kun-Rodrigues C, et al. **Investigating the genetic architecture of dementia with Lewy bodies: a two-stage genome-wide association study.** *Lancet Neurol* 2018. — APOE, SNCA and GBA1 are the three principal loci in DLB.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Investigating+the+genetic+architecture+of+dementia+with+Lewy+bodies+a+two-stage+genome-wide+association+study>
42. Chia R, Sabir MS, Bandres-Ciga S, et al. **Genome sequencing analysis identifies new loci associated with Lewy body dementia and provides insights into its genetic architecture.** *Nat Genet* 2021.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Genome+sequencing+analysis+identifies+new+loci+associated+with+Lewy+body+dementia+and+provides+insights+into+its+genetic+architecture>
43. Dickson DW, Heckman MG, Murray ME, et al. **APOE ε4 is associated with severity of Lewy body pathology independent of Alzheimer pathology.** *Neurology* 2018. — APOE4 worsens Lewy pathology independently of AD pathology. The model's dashed `APOE4` → `A_FIB` path.
    <https://pubmed.ncbi.nlm.nih.gov/?term=APOE+epsilon4+is+associated+with+severity+of+Lewy+body+pathology+independent+of+Alzheimer+pathology>
44. Singleton AB, Farrer M, Johnson J, et al. **α-Synuclein locus triplication causes Parkinson's disease.** *Science* 2003. — The gene dosage effect. `SNCADOSE`.
    <https://pubmed.ncbi.nlm.nih.gov/?term=alpha-Synuclein+locus+triplication+causes+Parkinson+disease>

## 7. The cholinergic system — presynaptic loss, postsynaptic preservation

45. Perry EK, Marshall E, Kerwin J, et al. **Evidence of a monoaminergic-cholinergic imbalance related to visual hallucinations in Lewy body dementia.** *J Neurochem* 1990. — Visual hallucinations and the cholinergic/monoaminergic imbalance.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Evidence+of+a+monoaminergic-cholinergic+imbalance+related+to+visual+hallucinations+in+Lewy+body+dementia>
46. Tiraboschi P, Hansen LA, Alford M, et al. **Cholinergic dysfunction in diseases with Lewy bodies.** *Neurology* 2000. — **The cortical ChAT decrement in DLB is larger than in AD.** The basis for the model's `NBM_0` (DLB 0.50 vs AD 0.72).
    <https://pubmed.ncbi.nlm.nih.gov/?term=Cholinergic+dysfunction+in+diseases+with+Lewy+bodies+Tiraboschi+Hansen+Alford>
47. Perry EK, Haroutunian V, Davis KL, et al. **Neocortical cholinergic activities differentiate Lewy body dementia from classical Alzheimer's disease.** *Neuroreport* 1994.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Neocortical+cholinergic+activities+differentiate+Lewy+body+dementia+from+classical+Alzheimer+disease>
48. Perry EK, Marshall E, Perry RH, et al. **Cholinergic and dopaminergic activities in senile dementia of Lewy body type.** *Alzheimer Dis Assoc Disord* 1990.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Cholinergic+and+dopaminergic+activities+in+senile+dementia+of+Lewy+body+type>
49. Shiozaki K, Iseki E, Uchiyama H, et al. **Alterations of muscarinic acetylcholine receptor subtypes in diffuse Lewy body disease: relation to Alzheimer's disease.** *J Neurol Neurosurg Psychiatry* 1999. — **Postsynaptic M1 is relatively preserved in DLB** — the direct basis for making `M1TGT = M1BASE − KTAUM1·PTAU` depend on tau alone in the model. The reason a ChEI works better in DLB than in AD.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Alterations+of+muscarinic+acetylcholine+receptor+subtypes+in+diffuse+Lewy+body+disease>
50. Ballard C, Piggott M, Johnson M, et al. **Delusions associated with elevated muscarinic binding in dementia with Lewy bodies.** *Ann Neurol* 2000. — Muscarinic binding is in fact raised.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Delusions+associated+with+elevated+muscarinic+binding+in+dementia+with+Lewy+bodies>
51. Pimlott SL, Piggott M, Owens J, et al. **Nicotinic acetylcholine receptor distribution in Alzheimer's disease, dementia with Lewy bodies, Parkinson's disease, and vascular dementia.** *Neuropsychopharmacology* 2004. — Thalamic and cortical reduction of α4β2 nAChR. `C_NIC`.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Nicotinic+acetylcholine+receptor+distribution+in+Alzheimer+disease+dementia+with+Lewy+bodies+Parkinson+disease+and+vascular+dementia>
52. Mesulam MM, Geula C. **Butyrylcholinesterase-rich neurons of the human cerebral cortex.** *Ann Neurol* 1988. — The cortical distribution of BuChE. The model's `BCHEFR`/`KBCHUP`: as AChE is lost, BuChE comes to carry a larger share of the hydrolysis, which is the basis for the prediction that a dual inhibitor (rivastigmine) holds a relative advantage in the later stages.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Butyrylcholinesterase-rich+neurons+of+the+human+cerebral+cortex+Mesulam+Geula>
53. Greig NH, Utsuki T, Ingram DK, et al. **Selective butyrylcholinesterase inhibition elevates brain acetylcholine, augments learning and lowers Alzheimer β-amyloid peptide in rodent.** *PNAS* 2005.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Selective+butyrylcholinesterase+inhibition+elevates+brain+acetylcholine+augments+learning+and+lowers+Alzheimer+beta-amyloid+peptide>

## 8. Dopaminergic system · DAT imaging · neuroleptic sensitivity

54. McKeith I, Fairbairn A, Perry R, et al. **Neuroleptic sensitivity in patients with senile dementia of Lewy body type.** *BMJ* 1992. — **The original paper.** A severe neuroleptic sensitivity reaction occurs in about half, with mortality 2–3 times higher. The calibration target for the model's `NSENS`.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Neuroleptic+sensitivity+in+patients+with+senile+dementia+of+Lewy+body+type+McKeith+Fairbairn+Perry>
55. Aarsland D, Perry R, Larsen JP, et al. **Neuroleptic sensitivity in Parkinson's disease and parkinsonian dementias.** *J Clin Psychiatry* 2005. — **The same sensitivity appears in PDD as well.** The reason it is right for the model to predict DLB and PDD similarly (and starkly differently from AD).
    <https://pubmed.ncbi.nlm.nih.gov/?term=Neuroleptic+sensitivity+in+Parkinson+disease+and+parkinsonian+dementias+Aarsland+Perry+Larsen>
56. Piggott MA, Marshall EF, Thomas N, et al. **Striatal dopaminergic markers in dementia with Lewy bodies, Alzheimer's and Parkinson's diseases: rostrocaudal distribution.** *Brain* 1999. — **Striatal D2 receptors are not upregulated in DLB** (whereas in PD they are). The direct basis for the model's `UPCAP = exp(−KSUPP·OLIGL)` term, and this single term generates the neuroleptic sensitivity and the blunted levodopa response at the same time.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Striatal+dopaminergic+markers+in+dementia+with+Lewy+bodies+Alzheimer+and+Parkinson+diseases+rostrocaudal+distribution>
57. Piggott MA, Ballard CG, Rowan E, et al. **Selective loss of dopamine D2 receptors in temporal cortex in dementia with Lewy bodies, association with cognitive decline.** *Synapse* 2007.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Selective+loss+of+dopamine+D2+receptors+in+temporal+cortex+in+dementia+with+Lewy+bodies>
58. McKeith I, O'Brien J, Walker Z, et al. **Sensitivity and specificity of dopamine transporter imaging with ¹²³I-FP-CIT SPECT in dementia with Lewy bodies: a phase III, multicentre study.** *Lancet Neurol* 2007. — The diagnostic performance of DaTSCAN. The model's `DATSBR` mapping and `TERMEXP` exponent.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Sensitivity+and+specificity+of+dopamine+transporter+imaging+with+123I-FP-CIT+SPECT+in+dementia+with+Lewy+bodies+phase+III>
59. Walker Z, Costa DC, Walker RWH, et al. **Differentiation of dementia with Lewy bodies from Alzheimer's disease using a dopaminergic presynaptic ligand.** *J Neurol Neurosurg Psychiatry* 2002.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Differentiation+of+dementia+with+Lewy+bodies+from+Alzheimer+disease+using+a+dopaminergic+presynaptic+ligand>
60. Kaasinen V, Vahlberg T. **Striatal dopamine in Parkinson disease: a meta-analysis of imaging studies.** *Ann Neurol* 2017. — Loss of striatal terminals precedes and exceeds loss of cell bodies. The model's `TERMEXP = 1.6`.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Striatal+dopamine+in+Parkinson+disease+a+meta-analysis+of+imaging+studies+Kaasinen+Vahlberg>
61. Farde L, Nordström AL, Wiesel FA, et al. **Positron emission tomographic analysis of central D1 and D2 dopamine receptor occupancy in patients treated with classical neuroleptics and clozapine.** *Arch Gen Psychiatry* 1992. — A D2 occupancy of ~80% is the threshold for extrapyramidal symptoms. The model's `EC50D2` and the position of the motor-reserve sigmoid.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Positron+emission+tomographic+analysis+of+central+D1+and+D2+dopamine+receptor+occupancy+in+patients+treated+with+classical+neuroleptics+and+clozapine>
62. Kapur S, Zipursky R, Jones C, et al. **Relationship between dopamine D2 occupancy, clinical response, and side effects: a double-blind PET study of first-episode schizophrenia.** *Am J Psychiatry* 2000.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Relationship+between+dopamine+D2+occupancy+clinical+response+and+side+effects+a+double-blind+PET+study+of+first-episode+schizophrenia>
63. Kegeles LS, Slifstein M, Frankle WG, et al. **Dose-occupancy study of striatal and extrastriatal dopamine D2 receptors by aripiprazole in schizophrenia with PET and [¹⁸F]fallypride.** *Neuropsychopharmacology* 2008.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Dose-occupancy+study+of+striatal+and+extrastriatal+dopamine+D2+receptors+by+aripiprazole+in+schizophrenia>

## 9. The serotonergic system and 5-HT2A

64. Ballard C, Piggott M, Johnson M, et al. **5-HT2A receptor binding in the ventral visual pathway and visual hallucinations in dementia with Lewy bodies.** *(and companion analyses)* — **5-HT2A binding in BA18/19 is increased in DLB with visual hallucinations.** The basis for the model's `H2TGT = 1 + KH2UP·denervation + KH2NEO·FIBN`, and the only one of the three transmitter systems in which the postsynaptic receptor is **up**regulated.
    <https://pubmed.ncbi.nlm.nih.gov/?term=5-HT2A+receptor+binding+ventral+visual+pathway+visual+hallucinations+dementia+with+Lewy+bodies>
65. Cheng AV, Ferrier IN, Morris CM, et al. **Cortical serotonin-S2 receptor binding in Lewy body dementia, Alzheimer's and Parkinson's diseases.** *J Neurol Sci* 1991.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Cortical+serotonin-S2+receptor+binding+in+Lewy+body+dementia+Alzheimer+and+Parkinson+diseases>
66. Sharp SI, Ballard CG, Chen CPL-H, Francis PT. **Aggressive behavior and neuroleptic medication are associated with increased number of α1-adrenoceptors in patients with Lewy body disease.** *Am J Geriatr Psychiatry* 2007.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Aggressive+behavior+and+neuroleptic+medication+are+associated+with+increased+number+of+alpha1-adrenoceptors+in+patients+with+Lewy+body+disease>
67. Vanover KE, Davis RE. **Role of 5-HT2A receptor antagonists in the treatment of insomnia.** *Nat Sci Sleep* 2010. — The pharmacology of 5-HT2A inverse agonists (the concept of constitutive activity). The model's `HTCONST` (agonist-independent Gq tone) term.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Role+of+5-HT2A+receptor+antagonists+in+the+treatment+of+insomnia+Vanover+Davis>
68. Weiner DM, Burstein ES, Nash N, et al. **5-Hydroxytryptamine2A receptor inverse agonists as antipsychotics.** *J Pharmacol Exp Ther* 2001. — The pharmacological basis of the pimavanserin class.
    <https://pubmed.ncbi.nlm.nih.gov/?term=5-Hydroxytryptamine2A+receptor+inverse+agonists+as+antipsychotics+Weiner+Burstein+Nash>

## 10. The noradrenergic system and the locus coeruleus

69. Del Tredici K, Braak H. **Dysfunction of the locus coeruleus–norepinephrine system and related circuitry in Parkinson's disease-related dementia.** *J Neurol Neurosurg Psychiatry* 2013. — The locus coeruleus is among the earliest targets. The model's `LC_0` (DLB 0.62 vs AD 0.80) and `KLLC`.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Dysfunction+of+the+locus+coeruleus-norepinephrine+system+and+related+circuitry+in+Parkinson+disease-related+dementia>
70. Szot P, White SS, Greenup JL, et al. **Compensatory changes in the noradrenergic nervous system in the locus ceruleus and hippocampus of postmortem subjects with Alzheimer's disease and dementia with Lewy bodies.** *J Neurosci* 2006.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Compensatory+changes+in+the+noradrenergic+nervous+system+in+the+locus+ceruleus+and+hippocampus+of+postmortem+subjects+with+Alzheimer+disease+and+dementia+with+Lewy+bodies>
71. Aston-Jones G, Cohen JD. **An integrative theory of locus coeruleus–norepinephrine function: adaptive gain and optimal performance.** *Annu Rev Neurosci* 2005. — The gain-control theory of LC-NE. The conceptual basis for the model's `NOISEG` (loss of LC+PPN raises the switching noise of the attentional state).
    <https://pubmed.ncbi.nlm.nih.gov/?term=An+integrative+theory+of+locus+coeruleus-norepinephrine+function+adaptive+gain+and+optimal+performance>

## 11. Cognitive fluctuation · EEG · thalamocortical circuitry

72. Bonanni L, Thomas A, Tiraboschi P, et al. **EEG comparisons in early Alzheimer's disease, dementia with Lewy bodies and Parkinson's disease with dementia patients with a 2-year follow-up.** *Brain* 2008. — The fall in dominant frequency in DLB and its **high temporal variability** (pre-alpha 5.6–7.9 Hz). The model's `EEGF` mapping.
    <https://pubmed.ncbi.nlm.nih.gov/?term=EEG+comparisons+in+early+Alzheimer+disease+dementia+with+Lewy+bodies+and+Parkinson+disease+with+dementia+patients+with+a+2-year+follow-up>
73. Bonanni L, Franciotti R, Nobili F, et al. **EEG markers of dementia with Lewy bodies: a multicenter cohort study.** *J Alzheimers Dis* 2016.
    <https://pubmed.ncbi.nlm.nih.gov/?term=EEG+markers+of+dementia+with+Lewy+bodies+a+multicenter+cohort+study>
74. Walker MP, Ayre GA, Cummings JL, et al. **Quantifying fluctuation in dementia with Lewy bodies, Alzheimer's disease, and vascular dementia.** *Neurology* 2000. — The quantitative comparison showing that fluctuation is characteristically large in DLB.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Quantifying+fluctuation+in+dementia+with+Lewy+bodies+Alzheimer+disease+and+vascular+dementia>
75. Ballard C, Walker M, O'Brien J, et al. **The characterisation and impact of 'fluctuating' cognition in dementia with Lewy bodies and Alzheimer's disease.** *Int J Geriatr Psychiatry* 2001.
    <https://pubmed.ncbi.nlm.nih.gov/?term=The+characterisation+and+impact+of+fluctuating+cognition+in+dementia+with+Lewy+bodies+and+Alzheimer+disease>
76. Delli Pizzi S, Franciotti R, Taylor JP, et al. **Thalamic involvement in fluctuating cognition in dementia with Lewy bodies: magnetic resonance evidences.** *Cereb Cortex* 2015. — Thalamic involvement. The model's `E_THAL` / `WTHAL`.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Thalamic+involvement+in+fluctuating+cognition+in+dementia+with+Lewy+bodies+magnetic+resonance+evidences>
77. Peraza LR, Cromarty R, Kobeleva X, et al. **Electroencephalographic derived network differences in Lewy body dementia compared to Alzheimer's disease patients.** *Sci Rep* 2018.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Electroencephalographic+derived+network+differences+in+Lewy+body+dementia+compared+to+Alzheimer+disease+patients>
78. Deco G, Jirsa VK, McIntosh AR. **Emerging concepts for the dynamical organization of resting-state activity in the brain.** *Nat Rev Neurosci* 2011. — Multistability of brain states and near-critical dynamics. The theoretical basis for the model's cubic attentional state (`dxdt_ATTM`) and for the structural claim that "the variance grows before the mean does".
    <https://pubmed.ncbi.nlm.nih.gov/?term=Emerging+concepts+for+the+dynamical+organization+of+resting-state+activity+in+the+brain>
79. Scheffer M, Bascompte J, Brock WA, et al. **Early-warning signals for critical transitions.** *Nature* 2009. — The rise in variance and autocorrelation on approach to a critical transition (critical slowing down). The form of `CURV`, `BISTAB` and `FLUCTGT`.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Early-warning+signals+for+critical+transitions+Scheffer+Bascompte+Brock>

## 12. Visual hallucinations and the PAD model

80. Collerton D, Perry E, McKeith I. **Why people see things that are not there: a novel Perception and Attention Deficit model for recurrent complex visual hallucinations.** *Behav Brain Sci* 2005. — **The original source of the model's multiplicative PAD structure (`(1−BOTTOMUP)·(1−TOPDOWN)·5HT2A`).** The claim that a perceptual deficit and an attentional deficit are **both** required.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Why+people+see+things+that+are+not+there+a+novel+Perception+and+Attention+Deficit+model+for+recurrent+complex+visual+hallucinations>
81. Shine JM, Halliday GM, Naismith SL, Lewis SJG. **Visual misperceptions and hallucinations in Parkinson's disease: dysfunction of attentional control networks?** *Mov Disord* 2011.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Visual+misperceptions+and+hallucinations+in+Parkinson+disease+dysfunction+of+attentional+control+networks>
82. Onofrj M, Taylor JP, Monaco D, et al. **Visual hallucinations in PD and Lewy body dementias: old and new hypotheses.** *Behav Neurol* 2013.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Visual+hallucinations+in+PD+and+Lewy+body+dementias+old+and+new+hypotheses>
83. Taylor JP, Firbank MJ, He J, et al. **Visual cortex in dementia with Lewy bodies: magnetic resonance imaging study.** *Br J Psychiatry* 2012.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Visual+cortex+in+dementia+with+Lewy+bodies+magnetic+resonance+imaging+study>
84. Lobotesis K, Fenwick JD, Phipps A, et al. **Occipital hypoperfusion on SPECT in dementia with Lewy bodies but not AD.** *Neurology* 2001. — Occipital hypometabolism/hypoperfusion. The model's `BOTTOMUP` decrement term.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Occipital+hypoperfusion+on+SPECT+in+dementia+with+Lewy+bodies+but+not+AD>
85. Lim SM, Katsifis A, Villemagne VL, et al. **The ¹⁸F-FDG PET cingulate island sign and comparison to ¹²³I-β-CIT SPECT for diagnosis of dementia with Lewy bodies.** *J Nucl Med* 2009. — cingulate island sign.
    <https://pubmed.ncbi.nlm.nih.gov/?term=The+18F-FDG+PET+cingulate+island+sign+and+comparison+to+123I-beta-CIT+SPECT+for+diagnosis+of+dementia+with+Lewy+bodies>
86. Murakami H, Shiraishi T, Umehara T, et al. **Retinal thinning and visual hallucination in dementia with Lewy bodies** *(OCT studies of inner retinal layers)*. — Thinning of the retinal layers. `V_RET`.
    <https://pubmed.ncbi.nlm.nih.gov/?term=retinal+thinning+optical+coherence+tomography+dementia+with+Lewy+bodies+visual+hallucination>

## 13. REM sleep behaviour disorder and sleep

87. Boeve BF, Silber MH, Ferman TJ, et al. **Clinicopathologic correlations in 172 cases of rapid eye movement sleep behavior disorder with or without a coexisting neurologic disorder.** *Sleep Med* 2013.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Clinicopathologic+correlations+in+172+cases+of+rapid+eye+movement+sleep+behavior+disorder+with+or+without+a+coexisting+neurologic+disorder>
88. Postuma RB, Iranzo A, Hu M, et al. **Risk and predictors of dementia and parkinsonism in idiopathic REM sleep behaviour disorder: a multicentre study.** *Brain* 2019. — An annual phenoconversion rate of about 6%. The basis for the prodromal cohort of scenario 19.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Risk+and+predictors+of+dementia+and+parkinsonism+in+idiopathic+REM+sleep+behaviour+disorder+a+multicentre+study>
89. Luppi PH, Clément O, Sapin E, et al. **The neuronal network responsible for paradoxical sleep and its dysfunctions causing narcolepsy and rapid eye movement sleep behavior disorder.** *Sleep Med Rev* 2011. — The SLD and GiV glycine/GABA circuitry. `Z_SLD`, `Z_GIV`.
    <https://pubmed.ncbi.nlm.nih.gov/?term=The+neuronal+network+responsible+for+paradoxical+sleep+and+its+dysfunctions+causing+narcolepsy+and+rapid+eye+movement+sleep+behavior+disorder>
90. Kunz D, Mahlberg R. **A two-part, double-blind, placebo-controlled trial of exogenous melatonin in REM sleep behaviour disorder.** *J Sleep Res* 2010. — `MELON`/`MELEFF`.
    <https://pubmed.ncbi.nlm.nih.gov/?term=A+two-part+double-blind+placebo-controlled+trial+of+exogenous+melatonin+in+REM+sleep+behaviour+disorder>
91. Ferman TJ, Smith GE, Dickson DW, et al. **Abnormal daytime sleepiness in dementia with Lewy bodies compared to Alzheimer's disease using the Multiple Sleep Latency Test.** *Alzheimers Res Ther* 2014. — Excessive daytime sleepiness. `Z_EDS`, `EDSS`.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Abnormal+daytime+sleepiness+in+dementia+with+Lewy+bodies+compared+to+Alzheimer+disease+using+the+Multiple+Sleep+Latency+Test>
92. Fronczek R, van Geest S, Frölich M, et al. **Hypocretin (orexin) loss in Alzheimer's disease.** *Neurobiol Aging* 2012. — Loss of orexin neurons. `Z_ORX`.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Hypocretin+orexin+loss+in+Alzheimer+disease+Fronczek>

## 14. Autonomic function and MIBG

93. Yoshita M, Arai H, Arai H, et al. **Diagnostic accuracy of ¹²³I-meta-iodobenzylguanidine myocardial scintigraphy in dementia with Lewy bodies: a multicenter study.** *PLoS One* 2015. — Sensitivity ~69%, specificity ~87%. The model's `MIBG = 1.05 + 1.75·CSYM` and `CSYM_0` (DLB 0.35 vs AD 0.85).
    <https://pubmed.ncbi.nlm.nih.gov/?term=Diagnostic+accuracy+of+123I-meta-iodobenzylguanidine+myocardial+scintigraphy+in+dementia+with+Lewy+bodies+a+multicenter+study>
94. Orimo S, Uchihara T, Nakamura A, et al. **Axonal α-synuclein aggregates herald centripetal degeneration of cardiac sympathetic nerve in Parkinson's disease.** *Brain* 2008. — α-Synuclein in the cardiac sympathetic terminals.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Axonal+alpha-synuclein+aggregates+herald+centripetal+degeneration+of+cardiac+sympathetic+nerve+in+Parkinson+disease>
95. Kaufmann H, Norcliffe-Kaufmann L, Palma JA, et al. **Natural history of pure autonomic failure: a United States prospective cohort.** *Ann Neurol* 2017.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Natural+history+of+pure+autonomic+failure+a+United+States+prospective+cohort>
96. Kaufmann H, Freeman R, Biaggioni I, et al. **Droxidopa for neurogenic orthostatic hypotension: a randomized, placebo-controlled, phase 3 trial.** *Neurology* 2014. — `DROXON`/`DROXEFF`, and the price paid in supine hypertension (`DROXSUP`).
    <https://pubmed.ncbi.nlm.nih.gov/?term=Droxidopa+for+neurogenic+orthostatic+hypotension+a+randomized+placebo-controlled+phase+3+trial>
97. Robertson AD, Udow SJ, Espay AJ, et al. **Orthostatic hypotension and dementia incidence: links and implications.** *Neuropsychiatr Dis Treat* 2019. — The mechanism (cerebral perfusion) by which OH is linked to cognitive decline. The model's `WOH` term.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Orthostatic+hypotension+and+dementia+incidence+links+and+implications>

## 15. Alzheimer co-pathology

98. Irwin DJ, Lee VMY, Trojanowski JQ. **Parkinson's disease dementia: convergence of α-synuclein, tau and amyloid-β pathologies.** *Nat Rev Neurosci* 2013. — The interaction of the three pathologies. The model's bilinear `B_CROSS` and `XSEED` terms.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Parkinson+disease+dementia+convergence+of+alpha-synuclein+tau+and+amyloid-beta+pathologies>
99. Masliah E, Rockenstein E, Veinbergs I, et al. **β-Amyloid peptides enhance α-synuclein accumulation and neuronal deficits in a transgenic mouse model linking Alzheimer's disease and Parkinson's disease.** *PNAS* 2001. — Direct experimental evidence for cross-seeding.
    <https://pubmed.ncbi.nlm.nih.gov/?term=beta-Amyloid+peptides+enhance+alpha-synuclein+accumulation+and+neuronal+deficits+in+a+transgenic+mouse+model+linking+Alzheimer+disease+and+Parkinson+disease>
100. Ferman TJ, Aoki N, Crook JE, et al. **The limbic and neocortical contribution of α-synuclein, tau, and amyloid-β to disease duration and severity in dementia with Lewy bodies.** *Alzheimers Dement* 2018. — Co-pathology changes both the phenotype and the disease duration. `B_PHENO`.
     <https://pubmed.ncbi.nlm.nih.gov/?term=The+limbic+and+neocortical+contribution+of+alpha-synuclein+tau+and+amyloid-beta+to+disease+duration+and+severity+in+dementia+with+Lewy+bodies>
101. Lemstra AW, de Beer MH, Teunissen CE, et al. **Concomitant AD pathology affects clinical manifestation and survival in dementia with Lewy bodies.** *J Neurol Neurosurg Psychiatry* 2017. — Fewer visual hallucinations and shorter survival. The model's dashed `B_PHENO → F_SURV` and `B_PHENO ⊣ V_VH` paths.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Concomitant+AD+pathology+affects+clinical+manifestation+and+survival+in+dementia+with+Lewy+bodies>

## 16. Neuroinflammation

102. Kim C, Ho DH, Suk JE, et al. **Neuron-released oligomeric α-synuclein is an endogenous agonist of TLR2 for paracrine activation of microglia.** *Nat Commun* 2013. — `I_MG`, `KMGON`.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Neuron-released+oligomeric+alpha-synuclein+is+an+endogenous+agonist+of+TLR2+for+paracrine+activation+of+microglia>
103. Liddelow SA, Guttenplan KA, Clarke LE, et al. **Neurotoxic reactive astrocytes are induced by activated microglia.** *Nature* 2017. — A1 astrocytes induced by C1q+IL-1α+TNF. `I_A1`.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Neurotoxic+reactive+astrocytes+are+induced+by+activated+microglia>
104. Sulzer D, Alcalay RN, Garretti F, et al. **T cells from patients with Parkinson's disease recognize α-synuclein peptides.** *Nature* 2017. — `I_TCELL`.
     <https://pubmed.ncbi.nlm.nih.gov/?term=T+cells+from+patients+with+Parkinson+disease+recognize+alpha-synuclein+peptides>
105. Surendranathan A, Su L, Mak E, et al. **Early microglial activation and peripheral inflammation in dementia with Lewy bodies.** *Brain* 2018. — **Activated early and exhausted late** — the sign reversal in the model's `PHAGO = MGA·(1 − KEXH·MGEXH)`.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Early+microglial+activation+and+peripheral+inflammation+in+dementia+with+Lewy+bodies>

## 17. Fluid and imaging biomarkers

106. Fairfoul G, McGuire LI, Pal S, et al. **α-Synuclein RT-QuIC in the CSF of patients with alpha-synucleinopathies.** *Ann Clin Transl Neurol* 2016. — The seed amplification assay. `W_SAA`.
     <https://pubmed.ncbi.nlm.nih.gov/?term=alpha-Synuclein+RT-QuIC+in+the+CSF+of+patients+with+alpha-synucleinopathies>
107. Siderowf A, Concha-Marambio L, Lafontant DE, et al. **Assessment of heterogeneity among participants in the Parkinson's Progression Markers Initiative cohort using α-synuclein seed amplification: a cross-sectional study.** *Lancet Neurol* 2023. — Sensitivity and specificity, and the differences between GBA/LRRK2 subtypes.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Assessment+of+heterogeneity+among+participants+in+the+Parkinson+Progression+Markers+Initiative+cohort+using+alpha-synuclein+seed+amplification>
108. Wang Z, Becker K, Donadio V, et al. **Skin α-synuclein aggregation seeding activity as a novel biomarker for Parkinson disease.** *JAMA Neurol* 2021.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Skin+alpha-synuclein+aggregation+seeding+activity+as+a+novel+biomarker+for+Parkinson+disease>
109. Pilotto A, Imarisio A, Carrarini C, et al. **Plasma neurofilament light chain predicts cognitive progression in prodromal and clinical dementia with Lewy bodies.** *J Alzheimers Dis* 2021. — `W_NFL`.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Plasma+neurofilament+light+chain+predicts+cognitive+progression+in+prodromal+and+clinical+dementia+with+Lewy+bodies>
110. Palmqvist S, Janelidze S, Quiroz YT, et al. **Discriminative accuracy of plasma phospho-tau217 for Alzheimer disease vs other neurodegenerative disorders.** *JAMA* 2020. — `B_PTAU`, and stratification by co-pathology.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Discriminative+accuracy+of+plasma+phospho-tau217+for+Alzheimer+disease+vs+other+neurodegenerative+disorders>

## 18. Cholinesterase inhibitor trials

111. McKeith I, Del Ser T, Spano P, et al. **Efficacy of rivastigmine in dementia with Lewy bodies: a randomised, double-blind, placebo-controlled international study.** *Lancet* 2000. — **The first large RCT in DLB.** Rivastigmine 6–12 mg/day, about 30% improvement in NPI-4 at 20 weeks; the effect is larger on the neuropsychiatric and fluctuation measures than on cognition. The calibration anchor for the model's prediction that "a ChEI moves the variance before it moves the mean".
     <https://pubmed.ncbi.nlm.nih.gov/?term=Efficacy+of+rivastigmine+in+dementia+with+Lewy+bodies+a+randomised+double-blind+placebo-controlled+international+study>
112. Mori E, Ikeda M, Kosaka K; Donepezil-DLB Study Investigators. **Donepezil for dementia with Lewy bodies: a randomized, placebo-controlled trial.** *Ann Neurol* 2012. — About +2.2 on MMSE with donepezil 10 mg, and improvement in NPI-2. The calibration anchor for the model's scenario 06.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Donepezil+for+dementia+with+Lewy+bodies+a+randomized+placebo-controlled+trial+Mori+Ikeda+Kosaka>
113. Ikeda M, Mori E, Matsuo K, et al. **Donepezil for dementia with Lewy bodies: a randomized, placebo-controlled, confirmatory phase III trial.** *Alzheimers Res Ther* 2015.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Donepezil+for+dementia+with+Lewy+bodies+a+randomized+placebo-controlled+confirmatory+phase+III+trial>
114. Stinton C, McKeith I, Taylor JP, et al. **Pharmacological management of Lewy body dementia: a systematic review and meta-analysis.** *Am J Psychiatry* 2015. — A synthesis of pharmacotherapy in DLB/PDD. The conclusion that ChEIs carry the most secure evidence.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Pharmacological+management+of+Lewy+body+dementia+a+systematic+review+and+meta-analysis>
115. Matsunaga S, Kishi T, Iwata N. **Combination therapy with cholinesterase inhibitors and memantine for Lewy body disorders: a systematic review and meta-analysis.** *Int J Geriatr Psychiatry* 2016. — Scenario 07.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Combination+therapy+with+cholinesterase+inhibitors+and+memantine+for+Lewy+body+disorders+a+systematic+review+and+meta-analysis>
116. Winblad B, Cummings J, Andreasen N, et al. **A six-month double-blind, randomized, placebo-controlled study of a transdermal patch in Alzheimer's disease — rivastigmine patch versus capsule (IDEAL).** *Int J Geriatr Psychiatry* 2007. — **The patch has efficacy similar to the capsule with fewer gastrointestinal adverse effects.** The result reproduced by the model's `FPATCH`/`GIAE` (Cmax-driven) structure.
     <https://pubmed.ncbi.nlm.nih.gov/?term=A+six-month+double-blind+randomized+placebo-controlled+study+of+a+transdermal+patch+in+Alzheimer+disease+rivastigmine+patch+versus+capsule+IDEAL>
117. Gobburu JVS, Tammara V, Lesko L, et al. **Pharmacokinetic-pharmacodynamic modeling of rivastigmine, a cholinesterase inhibitor, in patients with Alzheimer's disease.** *J Clin Pharmacol* 2001. — The separation of the plasma half-life (~1.5 h) from the **pharmacodynamic half-life (enzyme resynthesis, ~8–10 h).** The direct basis for the model holding carbamylated enzyme as separate states (`CARBA`/`CARBB`).
     <https://pubmed.ncbi.nlm.nih.gov/?term=Pharmacokinetic-pharmacodynamic+modeling+of+rivastigmine+a+cholinesterase+inhibitor+in+patients+with+Alzheimer+disease>
118. Bar-On P, Millard CB, Harel M, et al. **Kinetic and structural studies on the interaction of cholinesterases with the anti-Alzheimer drug rivastigmine.** *Biochemistry* 2002. — The pseudo-irreversible carbamylation mechanism.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Kinetic+and+structural+studies+on+the+interaction+of+cholinesterases+with+the+anti-Alzheimer+drug+rivastigmine>
119. Rogers SL, Cooper NM, Sukovaty R, et al. **Pharmacokinetic and pharmacodynamic profile of donepezil HCl following multiple oral doses.** *Br J Clin Pharmacol* 1998. — Half-life ~70 h, steady-state concentrations, and the fractional erythrocyte AChE inhibition. Calibration of the model's `CLDON`/`VDON`/`IC50DON`.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Pharmacokinetic+and+pharmacodynamic+profile+of+donepezil+HCl+following+multiple+oral+doses>
120. Kaasinen V, Någren K, Järvenpää T, et al. **Regional effects of donepezil and rivastigmine on cortical acetylcholinesterase activity in Alzheimer's disease.** *J Clin Psychopharmacol* 2002. — The in vivo fractional cortical AChE inhibition (of the order of 30–40%). Validation of the model's `INHACHE` range.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Regional+effects+of+donepezil+and+rivastigmine+on+cortical+acetylcholinesterase+activity+in+Alzheimer+disease>

## 19. Memantine

121. Emre M, Tsolaki M, Bonuccelli U, et al. **Memantine for patients with Parkinson's disease dementia or dementia with Lewy bodies: a randomised, double-blind, placebo-controlled trial.** *Lancet Neurol* 2010. — A small benefit on global impression. Scenario 07.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Memantine+for+patients+with+Parkinson+disease+dementia+or+dementia+with+Lewy+bodies+a+randomised+double-blind+placebo-controlled+trial>
122. Aarsland D, Ballard C, Walker Z, et al. **Memantine in patients with Parkinson's disease dementia or dementia with Lewy bodies: a double-blind, placebo-controlled, multicentre trial.** *Lancet Neurol* 2009.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Memantine+in+patients+with+Parkinson+disease+dementia+or+dementia+with+Lewy+bodies+a+double-blind+placebo-controlled+multicentre+trial>
123. Parsons CG, Stöffler A, Danysz W. **Memantine: a NMDA receptor antagonist that improves memory by restoration of homeostasis in the glutamatergic system — too little activation is bad, too much is even worse.** *Neuropharmacology* 2007. — Non-competitive blockade with a fast off-rate. `IC50MEM`, `EMAXMEM`.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Memantine+a+NMDA+receptor+antagonist+that+improves+memory+by+restoration+of+homeostasis+in+the+glutamatergic+system>

## 20. Antipsychotics and pimavanserin

124. Cummings J, Isaacson S, Mills R, et al. **Pimavanserin for patients with Parkinson's disease psychosis: a randomised, placebo-controlled phase 3 trial.** *Lancet* 2014. — About −3.06 on SAPS-PD. The calibration anchor for the model's scenario 10 (hallucination burden about −37%).
     <https://pubmed.ncbi.nlm.nih.gov/?term=Pimavanserin+for+patients+with+Parkinson+disease+psychosis+a+randomised+placebo-controlled+phase+3+trial>
125. Tariot PN, Cummings JL, Soto-Martin ME, et al. **Trial of pimavanserin in dementia-related psychosis (HARMONY).** *N Engl J Med* 2021. — A relapse-prevention design, including a DLB subgroup. Scenario 10.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Trial+of+pimavanserin+in+dementia-related+psychosis+HARMONY+Tariot+Cummings>
126. Nasrallah HA, Fedora R, Morton R. **Successful treatment of clozapine-nonresponsive refractory hallucinations and delusions with pimavanserin, a serotonin 5HT-2A receptor inverse agonist.** *Schizophr Res* 2019.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Successful+treatment+of+clozapine-nonresponsive+refractory+hallucinations+and+delusions+with+pimavanserin>
127. Ballard C, Kreitzman DL, Isaacson S, et al. **Long-term evaluation of open-label pimavanserin safety and tolerability in Parkinson's disease psychosis.** *Parkinsonism Relat Disord* 2020. — The magnitude of QTc prolongation. `QTSLOPE`.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Long-term+evaluation+of+open-label+pimavanserin+safety+and+tolerability+in+Parkinson+disease+psychosis>
128. Kyle K, Bronstein JM. **Treatment of psychosis in Parkinson's disease and dementia with Lewy bodies: a review.** *Parkinsonism Relat Disord* 2020.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Treatment+of+psychosis+in+Parkinson+disease+and+dementia+with+Lewy+bodies+a+review+Kyle+Bronstein>
129. Kurlan R, Cummings J, Raman R, Thal L; Alzheimer's Disease Cooperative Study Group. **Quetiapine for agitation or psychosis in patients with dementia and parkinsonism.** *Neurology* 2007. — The limited efficacy of quetiapine. Scenario 11.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Quetiapine+for+agitation+or+psychosis+in+patients+with+dementia+and+parkinsonism+Kurlan+Cummings>
130. Schneider LS, Dagerman KS, Insel P. **Risk of death with atypical antipsychotic drug treatment for dementia: meta-analysis of randomized placebo-controlled trials.** *JAMA* 2005. — The raised risk of death (the basis for the boxed warning). The model's multiplicative `BNS` hazard term.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Risk+of+death+with+atypical+antipsychotic+drug+treatment+for+dementia+meta-analysis+of+randomized+placebo-controlled+trials>
131. Ballard C, Hanney ML, Theodoulou M, et al. **The dementia antipsychotic withdrawal trial (DART-AD): long-term follow-up of a randomised placebo-controlled trial.** *Lancet Neurol* 2009.
     <https://pubmed.ncbi.nlm.nih.gov/?term=The+dementia+antipsychotic+withdrawal+trial+DART-AD+long-term+follow-up+of+a+randomised+placebo-controlled+trial>

## 21. Motor symptomatic therapy

132. Molloy S, McKeith IG, O'Brien JT, Burn DJ. **The role of levodopa in the management of dementia with Lewy bodies.** *J Neurol Neurosurg Psychiatry* 2005. — **The levodopa response is distinctly blunted in DLB** (typically a meaningful response in only about one third). The comparator for scenarios 15/16 and for `responder_rates()`.
     <https://pubmed.ncbi.nlm.nih.gov/?term=The+role+of+levodopa+in+the+management+of+dementia+with+Lewy+bodies>
133. Goetz CG, Blasucci LM, Leurgans S, Pappert EJ. **Olanzapine and clozapine: comparative effects on motor function in hallucinating PD patients.** *Neurology* 2000.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Olanzapine+and+clozapine+comparative+effects+on+motor+function+in+hallucinating+PD+patients>
134. Murata M, Odawara T, Hasegawa K, et al. **Adjunct zonisamide to levodopa for DLB parkinsonism: a randomized double-blind phase 2 study.** *Neurology* 2018. — Zonisamide, approved in Japan for the parkinsonism of DLB. Scenario 17.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Adjunct+zonisamide+to+levodopa+for+DLB+parkinsonism+a+randomized+double-blind+phase+2+study>
135. Murata M, Odawara T, Hasegawa K, et al. **Effect of zonisamide on parkinsonism in patients with dementia with Lewy bodies: a randomized, double-blind, placebo-controlled trial.** *Parkinsonism Relat Disord* 2020.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Effect+of+zonisamide+on+parkinsonism+in+patients+with+dementia+with+Lewy+bodies+a+randomized+double-blind+placebo-controlled+trial>
136. Nutt JG, Woodward WR, Hammerstad JP, et al. **The 'on-off' phenomenon in Parkinson's disease: relation to levodopa absorption and transport.** *N Engl J Med* 1984. — LNAA competition and gastric emptying. The model's `U_ENS ⊣ X_LDC`.
     <https://pubmed.ncbi.nlm.nih.gov/?term=The+on-off+phenomenon+in+Parkinson+disease+relation+to+levodopa+absorption+and+transport>
137. Contin M, Martinelli P. **Pharmacokinetics of levodopa.** *J Neurol* 2010. — Half-life ~1.5 h (with a DDCI), bioavailability, volume of distribution. `KALD`/`FLD`/`VLD`/`CLLD`.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Pharmacokinetics+of+levodopa+Contin+Martinelli>
138. Nagai M, Kubota M, Kimura Y, et al. **Anticholinergic burden and cognition in older adults with dementia** *(anticholinergic burden reviews)*. — The basis for scenario 09 (`ANTICH`).
     <https://pubmed.ncbi.nlm.nih.gov/?term=anticholinergic+burden+cognitive+decline+dementia+older+adults+systematic+review>

## 22. Disease-modifying candidates

139. Pagano G, Taylor KI, Anzures-Cabrera J, et al. **Trial of prasinezumab in early-stage Parkinson's disease (PASADENA).** *N Engl J Med* 2022. — The PK of an anti-α-synuclein antibody (4500 mg IV q4w) and the failure of the primary endpoint. The basis for scenario 21, and the reason the model was made to predict a **small** disease-modifying effect.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Trial+of+prasinezumab+in+early-stage+Parkinson+disease+PASADENA>
140. Lang AE, Siderowf AD, Macklin EA, et al. **Trial of cinpanemab in early Parkinson's disease (SPARK).** *N Engl J Med* 2022.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Trial+of+cinpanemab+in+early+Parkinson+disease+SPARK>
141. Cole TA, Zhao H, Collier TJ, et al. **α-Synuclein antisense oligonucleotides as a disease-modifying therapy for Parkinson's disease.** *JCI Insight* 2021. — `Y_ASO`.
     <https://pubmed.ncbi.nlm.nih.gov/?term=alpha-Synuclein+antisense+oligonucleotides+as+a+disease-modifying+therapy+for+Parkinson+disease>
142. Jennings D, Huntwork-Rodriguez S, Vissers MFJM, et al. **LRRK2 inhibition by BIIB122 in healthy participants and patients with Parkinson's disease.** *Mov Disord* 2023. — `Y_LRRK`.
     <https://pubmed.ncbi.nlm.nih.gov/?term=LRRK2+inhibition+by+BIIB122+in+healthy+participants+and+patients+with+Parkinson+disease>
143. van Dyck CH, Swanson CJ, Aisen P, et al. **Lecanemab in early Alzheimer's disease.** *N Engl J Med* 2023. — A treatment acting only on the co-pathology arm (`Y_ADUC`), and the ARIA risk in APOE ε4.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Lecanemab+in+early+Alzheimer+disease+van+Dyck+Swanson+Aisen>

## 23. Prognosis and survival

144. Mueller C, Soysal P, Rongve A, et al. **Survival time and differences between dementia with Lewy bodies and Alzheimer's disease following diagnosis: a meta-analysis of longitudinal studies.** *Ageing Res Rev* 2019. — Survival in DLB is shorter than in AD. Calibration of the model's `H0`/`BCOG`/`BMOT`/`BFALL` (median survival about 4.4 years).
     <https://pubmed.ncbi.nlm.nih.gov/?term=Survival+time+and+differences+between+dementia+with+Lewy+bodies+and+Alzheimer+disease+following+diagnosis+a+meta-analysis+of+longitudinal+studies>
145. Williams MM, Xiong C, Morris JC, Galvin JE. **Survival and mortality differences between dementia with Lewy bodies vs Alzheimer disease.** *Neurology* 2006.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Survival+and+mortality+differences+between+dementia+with+Lewy+bodies+vs+Alzheimer+disease>
146. Rongve A, Soennesyn H, Skogseth R, et al. **Cognitive decline in dementia with Lewy bodies: a 5-year prospective cohort study.** *BMJ Open* 2016. — The annual rate of MMSE decline. The principal anchor for the model's natural-history calibration.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Cognitive+decline+in+dementia+with+Lewy+bodies+a+5-year+prospective+cohort+study+Rongve>
147. Kramberger MG, Auestad B, Garcia-Ptacek S, et al. **Long-term cognitive decline in dementia with Lewy bodies in a large multicenter, pooled cohort.** *J Alzheimers Dis* 2017.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Long-term+cognitive+decline+in+dementia+with+Lewy+bodies+in+a+large+multicenter+pooled+cohort>
148. Bostrom F, Jonsson L, Minthon L, Londos E. **Patients with dementia with Lewy bodies have more impaired quality of life than patients with Alzheimer disease.** *Alzheimer Dis Assoc Disord* 2007. — Carer burden. `F_CARE`.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Patients+with+dementia+with+Lewy+bodies+have+more+impaired+quality+of+life+than+patients+with+Alzheimer+disease>

## 24. QSP methodology

149. Baron KT, Elmokadem A, et al. **mrgsolve: Simulate from ODE-Based Models.** R package. — The simulation engine of this model.
     <https://mrgsolve.org/>
150. Gadkar K, Kirouac DC, Mager DE, et al. **A six-stage workflow for robust application of systems pharmacology.** *CPT Pharmacometrics Syst Pharmacol* 2016.
     <https://pubmed.ncbi.nlm.nih.gov/?term=A+six-stage+workflow+for+robust+application+of+systems+pharmacology>
151. Geerts H, Spiros A, Roberts P, Carr R. **Quantitative systems pharmacology as an extension of PK/PD modeling in CNS research and development.** *J Pharmacokinet Pharmacodyn* 2013. — The approach of representing receptor-level transmitters explicitly in CNS QSP. Directly adjacent to this model's three-transmitter structure.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Quantitative+systems+pharmacology+as+an+extension+of+PK+PD+modeling+in+CNS+research+and+development>
152. Geerts H, Roberts P, Spiros A, Carr R. **A strategy for developing new treatment paradigms for neuropsychiatric and neurocognitive symptoms in Alzheimer's disease.** *Front Pharmacol* 2013.
     <https://pubmed.ncbi.nlm.nih.gov/?term=A+strategy+for+developing+new+treatment+paradigms+for+neuropsychiatric+and+neurocognitive+symptoms+in+Alzheimer+disease>
153. Nicholas T, Duvvuri S, Leurent C, et al. **Systems pharmacology modeling in neuroscience: prediction and outcome of PF-04995274, a 5-HT4 partial agonist, in a clinical scopolamine impairment trial.** *Adv Alzheimers Dis* 2013. — A case of prospective validation of a QSP model in which the cholinergic transmitter is represented explicitly.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Systems+pharmacology+modeling+in+neuroscience+prediction+and+outcome+of+PF-04995274+a+5-HT4+partial+agonist+in+a+clinical+scopolamine+impairment+trial>
154. Bhattacharya S, Zhang Q, Andersen ME. **A deterministic map of Waddington's epigenetic landscape for cell fate decisions.** *BMC Syst Biol* 2011. — Biological modelling of the saddle-node bifurcation and of hysteresis. The model's GBA1 switch and cubic attentional state.
     <https://pubmed.ncbi.nlm.nih.gov/?term=A+deterministic+map+of+Waddington+epigenetic+landscape+for+cell+fate+decisions>

---

## What the literature does NOT support

Stated explicitly in the interest of honesty. The items below are **hypotheses of the model**, and at present have no evidential support.

1. **The specific functional form of `FLUCTGT`.** That cognitive fluctuation is the *variance* of a
   bistable attentional state is an **interpretation** consistent with the literature (#4, #74, #78, #79), but the
   functional form `SIGN0·(FLBASE + BISTAB·NOISEG)` has not been fitted to any data. The testable prediction
   this generates is that "a ChEI improves the fluctuation measure proportionally more than it improves MMSE",
   and the model computes that ratio as about 2.5-fold — the ratio itself is a falsifiable claim.
2. **`KSUPP` (the strength with which limbic α-synuclein suppresses postsynaptic striatal integrity).** #56
   supplies the observation that "D2 is not upregulated in DLB", but the **slope** of that relationship
   has never been measured. The model back-calculated this single parameter to match the incidence of
   neuroleptic sensitivity (#54, #55).
3. **The ambroxol "deadline".** The contrast between scenarios 19 and 20 rests entirely on the premise
   that the GBA1 feedback loop genuinely possesses a saddle-node bifurcation. #31 and #39 support the
   **existence** of the loop but do not demonstrate **bistability**. This is the largest extrapolation in
   this model, and at the same time its most easily testable prediction (it suffices to compare the same
   exposure in a prodromal against an established cohort).
4. **`MOTNDA` (the non-dopaminergic axial motor burden).** That a substantial part of the parkinsonism of
   DLB does not respond to levodopa (#132) is well established, but making that fraction a linear function
   of the neocortical and limbic fibril burden is the model's own choice.
5. **All of the neuronal loss rate constants (`KL*`).** These are not individually measured values but
   values tuned so as to match the natural-history calibration anchors (annual MMSE decline, annual
   MDS-UPDRS III rise, MIBG, DaTSCAN, median survival) simultaneously. Because they are correlated with one
   another, no biological meaning should be attached to the individual values.

---

*This reference list is for educational and research purposes. Do not use it for clinical decision-making.*
*This reference list is for education and research only. Not for clinical use.*
