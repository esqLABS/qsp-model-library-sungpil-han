# Congenital Hyperinsulinism (CHI) — References

The last line of each section states **which term of this model depends on that literature**.
The definitions of the numbers actually used in calibration (N1–N8, C1) and of the predictions (P1–P10) are in
the headers of [`README_en.md`](README_en.md) and [`chi_mrgsolve_model.R`](chi_mrgsolve_model.R).

---

## 1. Reviews · diagnosis · guidelines

1. Stanley CA. **Perspective on the genetics and diagnosis of congenital hyperinsulinism disorders.** J Clin Endocrinol Metab. 2016;101(3):815–26. <https://pubmed.ncbi.nlm.nih.gov/26908103/>
2. De Leon DD, Stanley CA. **Congenital hypoglycemia disorders: new aspects of etiology, diagnosis, treatment and outcomes.** Pediatr Diabetes. 2017;18(1):3–9. <https://pubmed.ncbi.nlm.nih.gov/28074627/>
3. Lord K, De León DD. **Monogenic hyperinsulinemic hypoglycemia: current insights into the pathogenesis and management.** Int J Pediatr Endocrinol. 2013;2013(1):3. <https://pubmed.ncbi.nlm.nih.gov/23448374/>
4. Thornton PS, Stanley CA, De Leon DD, et al. **Recommendations from the Pediatric Endocrine Society for evaluation and management of persistent hypoglycemia in neonates, infants, and children.** J Pediatr. 2015;167(2):238–45. <https://pubmed.ncbi.nlm.nih.gov/25957977/>
5. Banerjee I, Salomon-Estebanez M, Shah P, et al. **Therapies and outcomes of congenital hyperinsulinism-induced hypoglycaemia.** Diabet Med. 2019;36(1):9–21. <https://pubmed.ncbi.nlm.nih.gov/30246418/>
6. Demirbilek H, Hussain K. **Congenital hyperinsulinism: diagnosis and treatment update.** J Clin Res Pediatr Endocrinol. 2017;9(Suppl 2):69–87. <https://pubmed.ncbi.nlm.nih.gov/29280739/>
7. Rosenfeld E, Ganguly A, De Leon DD. **Congenital hyperinsulinism disorders: genetic and clinical characteristics.** Am J Med Genet C. 2019;181(4):682–92. <https://pubmed.ncbi.nlm.nih.gov/31414570/>
8. Galcheva S, Demirbilek H, Al-Khawaga S, Hussain K. **The genetic and molecular mechanisms of congenital hyperinsulinism.** Front Endocrinol. 2019;10:111. <https://pubmed.ncbi.nlm.nih.gov/30873120/>

> Model dependency: the working definition of CHI as GIR > 8 mg/kg/min (the diagnostic node), the critical sample items,
> and the management standard of a target glucose of 70 mg/dL — the model **does not assume** this target but
> re-derives it from the arithmetic of brain fuel (P5).

---

## 2. The K_ATP channel — the molecule that becomes the axis (g) of this model (ABCC8 / KCNJ11)

9. Thomas PM, Cote GJ, Wohllk N, et al. **Mutations in the sulfonylurea receptor gene in familial persistent hyperinsulinemic hypoglycemia of infancy.** Science. 1995;268(5209):426–9. <https://pubmed.ncbi.nlm.nih.gov/7716548/>
10. Nichols CG, Shyng SL, Nestorowicz A, et al. **Adenosine diphosphate as an intracellular regulator of insulin secretion.** Science. 1996;272(5269):1785–7. <https://pubmed.ncbi.nlm.nih.gov/8650576/>
11. Ashcroft FM. **ATP-sensitive potassium channelopathies: focus on insulin secretion.** J Clin Invest. 2005;115(8):2047–58. <https://pubmed.ncbi.nlm.nih.gov/16075048/>
12. Ashcroft FM, Rorsman P. **K_ATP channels and islet hormone secretion: new insights and controversies.** Nat Rev Endocrinol. 2013;9(11):660–9. <https://pubmed.ncbi.nlm.nih.gov/24042324/>
13. Shyng SL, Nichols CG. **Octameric stoichiometry of the K_ATP channel complex.** J Gen Physiol. 1997;110(6):655–64. <https://pubmed.ncbi.nlm.nih.gov/9382894/>
14. Cartier EA, Conti LR, Vandenberg CA, Shyng SL. **Defective trafficking and function of K_ATP channels caused by a sulfonylurea receptor 1 mutation associated with persistent hyperinsulinemic hypoglycemia of infancy.** Proc Natl Acad Sci USA. 2001;98(5):2882–7. <https://pubmed.ncbi.nlm.nih.gov/11226335/>
15. Snider KE, Becker S, Boyajian L, et al. **Genotype and phenotype correlations in 417 children with congenital hyperinsulinism.** J Clin Endocrinol Metab. 2013;98(2):E355–63. <https://pubmed.ncbi.nlm.nih.gov/23275527/>
16. Pinney SE, MacMullen C, Becker S, et al. **Clinical characteristics and biochemical mechanisms of congenital hyperinsulinism associated with dominant K_ATP channel mutations.** J Clin Invest. 2008;118(8):2877–86. <https://pubmed.ncbi.nlm.nih.gov/18596924/>
17. Macmullen CM, Zhou Q, Snider KE, et al. **Diazoxide-unresponsive congenital hyperinsulinism in children with dominant mutations of the beta-cell sulfonylurea receptor SUR1.** Diabetes. 2011;60(6):1797–804. <https://pubmed.ncbi.nlm.nih.gov/21536951/>
18. Saint-Martin C, Arnoux JB, de Lonlay P, Bellanné-Chantelot C. **KATP channel mutations in congenital hyperinsulinism.** Semin Pediatr Surg. 2011;20(1):18–22. <https://pubmed.ncbi.nlm.nih.gov/?term=KATP+channel+mutations+congenital+hyperinsulinism+Saint-Martin>

> Model dependency: **the core.** The structure in which the surface channel density g is the single lesion variable,
> and the point that a large number of recessive alleles drive g→0 through failure of *trafficking* (entry 14),
> is the molecular basis for "diazoxide can do nothing if there is no channel to open" (P2).
> Entries 15, 16, and 17 supply the genotype-phenotype separation, and the model **only reproduces** that separation;
> it does not calibrate to it.

---

## 3. Membrane potential · Ca²⁺ · secretion coupling (the demonstration of the voltage divider)

19. Ashcroft FM, Harrison DE, Ashcroft SJ. **Glucose induces closure of single potassium channels in isolated rat pancreatic beta-cells.** Nature. 1984;312(5993):446–8. <https://pubmed.ncbi.nlm.nih.gov/6095103/>
20. Rorsman P, Ashcroft FM. **Pancreatic β-cell electrical activity and insulin secretion: of mice and men.** Physiol Rev. 2018;98(1):117–214. <https://pubmed.ncbi.nlm.nih.gov/29212789/>
21. Henquin JC. **Triggering and amplifying pathways of regulation of insulin secretion by glucose.** Diabetes. 2000;49(11):1751–60. <https://pubmed.ncbi.nlm.nih.gov/11078440/>
22. Henquin JC. **Regulation of insulin secretion: a matter of phase control and amplitude modulation.** Diabetologia. 2009;52(5):739–51. <https://pubmed.ncbi.nlm.nih.gov/19288070/>
23. Straub SG, Sharp GW. **Hypothesis: one rate-limiting step controls the magnitude of both phases of glucose-stimulated insulin secretion.** Am J Physiol Cell Physiol. 2004;287(3):C565–71. <https://pubmed.ncbi.nlm.nih.gov/15308464/>
24. Matschinsky FM. **Banting Lecture 1995: A lesson in metabolic regulation inspired by the glucokinase glucose sensor paradigm.** Diabetes. 1996;45(2):223–41. <https://pubmed.ncbi.nlm.nih.gov/8549869/>

> Model dependency: the chain `Popen → G_KATP → V_m → fCa → secretion`,
> and the term in which the **amplifying pathway** (entries 21, 22) has a floor of 35 % in hypoglycaemia.
> That floor is the reason secretion at g=0 is finite rather than infinite.

---

## 4. Non-K_ATP genotypes — ways of touching the same channel from another side

25. Stanley CA, Lieu YK, Hsu BY, et al. **Hyperinsulinism and hyperammonemia in infants with regulatory mutations of the glutamate dehydrogenase gene.** N Engl J Med. 1998;338(19):1352–7. <https://pubmed.ncbi.nlm.nih.gov/9571255/>
26. Stanley CA, Fang J, Kutyna K, et al. **Molecular basis and characterization of the hyperinsulinism/hyperammonemia syndrome: predominance of mutations in exons 11 and 12 of the glutamate dehydrogenase gene.** Diabetes. 2000;49(4):667–73. <https://pubmed.ncbi.nlm.nih.gov/10871207/>
27. Kelly A, Ng D, Ferry RJ Jr, et al. **Acute insulin responses to leucine in children with the hyperinsulinism/hyperammonemia syndrome.** J Clin Endocrinol Metab. 2001;86(8):3724–8. <https://pubmed.ncbi.nlm.nih.gov/11502802/>
28. Glaser B, Kesavan P, Heyman M, et al. **Familial hyperinsulinism caused by an activating glucokinase mutation.** N Engl J Med. 1998;338(4):226–30. <https://pubmed.ncbi.nlm.nih.gov/9435328/>
29. Christesen HB, Jacobsen BB, Odili S, et al. **The second activating glucokinase mutation (A456V): implications for glucose homeostasis and diabetes therapy.** Diabetes. 2002;51(4):1240–6. <https://pubmed.ncbi.nlm.nih.gov/11916951/>
30. Clayton PT, Eaton S, Aynsley-Green A, et al. **Hyperinsulinism in short-chain L-3-hydroxyacyl-CoA dehydrogenase deficiency reveals the importance of beta-oxidation in insulin secretion.** J Clin Invest. 2001;108(3):457–65. <https://pubmed.ncbi.nlm.nih.gov/11489939/>
31. Pearson ER, Boj SF, Steele AM, et al. **Macrosomia and hyperinsulinaemic hypoglycaemia in patients with heterozygous mutations in the HNF4A gene.** PLoS Med. 2007;4(4):e118. <https://pubmed.ncbi.nlm.nih.gov/17407387/>
32. Otonkoski T, Jiao H, Kaminen-Ahola N, et al. **Physical exercise-induced hypoglycemia caused by failed silencing of monocarboxylate transporter 1 in pancreatic beta cells.** Am J Hum Genet. 2007;81(3):467–74. <https://pubmed.ncbi.nlm.nih.gov/17701894/>
33. González-Barroso MM, Giurgea I, Bouillaud F, et al. **Mutations in UCP2 in congenital hyperinsulinism reveal a role for regulation of insulin secretion.** PLoS One. 2008;3(12):e3850. <https://pubmed.ncbi.nlm.nih.gov/19065272/>
34. Flanagan SE, Vairo F, Johnson MB, et al. **A CACNA1D mutation in a patient with persistent hyperinsulinaemic hypoglycaemia, heart defects, and severe hypotonia.** Pediatr Diabetes. 2017;18(4):320–3. <https://pubmed.ncbi.nlm.nih.gov/28107785/>

> Model dependency: these genotypes **all keep g=1** and change only Rt (the metabolic signal), the S₀.₅ of GCK,
> or the GTP brake on GDH. The model therefore puts forward as a prediction that they
> **respond well** to diazoxide (P2, scenario 4).
> Entry 27 is the comparator for the leucine load prediction (P8, 37 vs 16 mg/dL).

---

## 5. Focal versus diffuse · pathology · ¹⁸F-DOPA PET

35. de Lonlay P, Fournet JC, Rahier J, et al. **Somatic deletion of the imprinted 11p15 region in sporadic persistent hyperinsulinemic hypoglycemia of infancy is specific of focal adenomatous hyperplasia.** J Clin Invest. 1997;100(4):802–7. <https://pubmed.ncbi.nlm.nih.gov/9259578/>
36. Rahier J, Guiot Y, Sempoux C. **Persistent hyperinsulinaemic hypoglycaemia of infancy: a heterogeneous syndrome unrelated to nesidioblastosis.** Arch Dis Child Fetal Neonatal Ed. 2000;82(2):F108–12. <https://pubmed.ncbi.nlm.nih.gov/10685983/>
37. Ribeiro MJ, De Lonlay P, Delzescaux T, et al. **Characterization of hyperinsulinism in infancy assessed with PET and ¹⁸F-fluoro-L-DOPA.** J Nucl Med. 2005;46(4):560–6. <https://pubmed.ncbi.nlm.nih.gov/15809475/>
38. Hardy OT, Hernandez-Pampaloni M, Saffer JR, et al. **Accuracy of ¹⁸F-fluorodopa positron emission tomography for diagnosing and localizing focal congenital hyperinsulinism.** J Clin Endocrinol Metab. 2007;92(12):4706–11. <https://pubmed.ncbi.nlm.nih.gov/17895314/>
39. Laje P, States LJ, Zhuang H, et al. **Accuracy of PET/CT scan in the diagnosis of the focal form of congenital hyperinsulinism.** J Pediatr Surg. 2013;48(2):388–93. <https://pubmed.ncbi.nlm.nih.gov/23414871/>
40. Sempoux C, Capito C, Bellanné-Chantelot C, et al. **Morphological mosaicism of the pancreatic islets: a novel anatomopathological form of persistent hyperinsulinemic hypoglycemia of infancy.** J Clin Endocrinol Metab. 2011;96(12):3785–93. <https://pubmed.ncbi.nlm.nih.gov/21956412/>

> Model dependency: the structure of two β-cell populations (`w_ab`, `dens_a`).
> The description in entries 36 and 40 of "high-density β cells with enlarged nuclei" meshes with the quantitative requirement of P10
> (lesion fraction × focal density ≈ 1.0). The model puts this forward as **both a prediction and a requirement**:
> for a small lesion to produce the same severity as diffuse disease, the density has to be ~10-fold.

---

## 6. Diazoxide — the drug that multiplies g

41. Drash A, Wolff F. **Drug therapy in leucine-sensitive hypoglycemia.** Metabolism. 1964;13:487–92. <https://pubmed.ncbi.nlm.nih.gov/14169215/>
42. Trube G, Rorsman P, Ohno-Shosaku T. **Opposite effects of tolbutamide and diazoxide on the ATP-dependent K⁺ channel in mouse pancreatic beta-cells.** Pflugers Arch. 1986;407(5):493–9. <https://pubmed.ncbi.nlm.nih.gov/3786110/>
43. Shyng SL, Ferrigni T, Nichols CG. **Regulation of K_ATP channel activity by diazoxide and MgADP: distinct functions of the two nucleotide binding folds of the sulfonylurea receptor.** J Gen Physiol. 1997;110(6):643–54. <https://pubmed.ncbi.nlm.nih.gov/9382893/>
44. Herman GA, Hussain K. **Diazoxide in the management of congenital hyperinsulinism: dose, response and adverse effects.** (a review of clinical series) Arch Dis Child. 2015;100(9):884–5. <https://pubmed.ncbi.nlm.nih.gov/25977564/>
45. Welters A, Lerch C, Kummer S, et al. **Long-term medical treatment in congenital hyperinsulinism: a descriptive analysis in a large cohort of patients from different clinical centers.** Orphanet J Rare Dis. 2015;10:150. <https://pubmed.ncbi.nlm.nih.gov/26608306/>
46. Thornton P, Truong L, Reynolds C, et al. **Rate of serious adverse events associated with diazoxide treatment of patients with hyperinsulinism.** Horm Res Paediatr. 2019;91(1):25–32. <https://pubmed.ncbi.nlm.nih.gov/30836359/>
47. Timlin MR, Black AB, Delaney HM, et al. **Development of pulmonary hypertension during treatment with diazoxide: a case series and literature review.** Pediatr Cardiol. 2017;38(6):1247–50. <https://pubmed.ncbi.nlm.nih.gov/28642988/>
48. Pruitt LG, Bloomstone J, et al. **Diazoxide pharmacokinetics in the newborn** — data on a neonatal t½ of 9–30 h and protein binding >90 %. Clin Pharmacol Ther. (a classic source) <https://pubmed.ncbi.nlm.nih.gov/?term=diazoxide+pharmacokinetics+newborn>

> Model dependency: **entry 43 is decisive.** Because diazoxide acts MgADP-dependently at NBD2 of SUR1,
> the model writes the drug not as "forcing a fixed proportion of channels open" but as
> **pushing the ATP/ADP set point of the channel (R50) to the right**. Written the first way to begin with,
> 2.5 mg/kg/day abolished secretion in the normal β cell (a glucose of 260 mg/dL); this was found and corrected during
> verification. The PK is entry 48, and the toxicity nodes are entries 46 and 47.

---

## 7. Somatostatin analogues — the drug that adds to the divider

49. Yamada Y, Post SR, Wang K, et al. **Cloning and functional characterization of a family of human and mouse somatostatin receptors expressed in brain, gastrointestinal tract, and kidney.** Proc Natl Acad Sci USA. 1992;89(1):251–5. <https://pubmed.ncbi.nlm.nih.gov/1346068/>
50. Rorsman P, Braun M. **Regulation of insulin secretion in human pancreatic islets.** Annu Rev Physiol. 2013;75:155–79. <https://pubmed.ncbi.nlm.nih.gov/22974438/>
51. Kailey B, van de Bunt M, Cheley S, et al. **SSTR2 is the functionally dominant somatostatin receptor in human pancreatic β- and α-cells.** Am J Physiol Endocrinol Metab. 2012;303(9):E1107–16. <https://pubmed.ncbi.nlm.nih.gov/22932785/>
52. Le Quan Sang KH, Arnoux JB, Mamoune A, et al. **Successful treatment of congenital hyperinsulinism with long-acting release octreotide.** Eur J Endocrinol. 2012;166(2):333–9. <https://pubmed.ncbi.nlm.nih.gov/22084154/>
53. Demirbilek H, Shah P, Arya VB, et al. **Long-term follow-up of children with congenital hyperinsulinism on octreotide therapy.** J Clin Endocrinol Metab. 2014;99(10):3660–7. <https://pubmed.ncbi.nlm.nih.gov/25004250/>
54. Laje P, Halaby L, Adzick NS, Stanley CA. **Necrotizing enterocolitis in neonates receiving octreotide for the management of congenital hyperinsulinism.** Pediatr Diabetes. 2010;11(2):142–7. <https://pubmed.ncbi.nlm.nih.gov/19558634/>
55. Corda H, Kummer S, Welters A, et al. **Treatment with long-acting lanreotide autogel in early infancy in patients with severe neonatal hyperinsulinism.** Orphanet J Rare Dis. 2017;12(1):108. <https://pubmed.ncbi.nlm.nih.gov/28615019/>

> Model dependency: **entry 51 is the second axis of this model.** Because SSTR2 lowers cAMP through Gi while
> at the same time **opening a G-protein-gated K⁺ conductance (GIRK)**, octreotide brings in
> a **second channel** that does not carry the mutation. That is why an effect remains even at g=0
> (P3: 59 %), and it is the only reason the conclusion comes out exactly opposite to diazoxide's.
> `gGIRK = 1.2` is the **one number** this model took from CHI (C1), and it was fitted to entry 53.
> Entry 54 is the neonatal NEC risk node.

---

## 8. Glucagon · insulin receptor targets · other agents

56. Neylon OM, Moran MM, Pellicano A, et al. **Successful subcutaneous glucagon use for persistent hypoglycaemia in congenital hyperinsulinism.** J Pediatr Endocrinol Metab. 2013;26(11-12):1157–61. <https://pubmed.ncbi.nlm.nih.gov/23751387/>
57. Hawkes CP, Lado JJ, Givler S, De Leon DD. **The effect of continuous intravenous glucagon on glucose requirements in infants with congenital hyperinsulinism.** JIMD Rep. 2019;45:45–50. <https://pubmed.ncbi.nlm.nih.gov/30569021/>
58. Ferrara CT, Boodhansingh KE, Paradies E, et al. **Novel hypoglycemia phenotype in congenital hyperinsulinism due to dominant mutations of uncoupling protein 2.** J Clin Endocrinol Metab. 2017;102(3):942–9. <https://pubmed.ncbi.nlm.nih.gov/28324055/>
59. Calabria AC, Li C, Gallagher PR, et al. **GLP-1 receptor antagonist exendin-(9-39) elevates fasting blood glucose levels in congenital hyperinsulinism owing to inactivating mutations in the ATP-sensitive K⁺ channel.** Diabetes. 2012;61(10):2585–91. <https://pubmed.ncbi.nlm.nih.gov/22855730/>
60. Ng CM, Tang F, Seeholzer SH, et al. **Population pharmacokinetics of exendin-(9-39) and clinical dose selection in patients with congenital hyperinsulinism.** Br J Clin Pharmacol. 2018;84(3):520–32. <https://pubmed.ncbi.nlm.nih.gov/29148594/>
61. Senniappan S, Alexandrescu S, Tatevian N, et al. **Sirolimus therapy in infants with severe hyperinsulinemic hypoglycemia.** N Engl J Med. 2014;370(12):1131–7. <https://pubmed.ncbi.nlm.nih.gov/24645944/>
62. Szymanowski M, Estebanez MS, Padidela R, et al. **mTOR inhibitors for the treatment of severe congenital hyperinsulinism: perspectives on limited therapeutic success.** J Clin Endocrinol Metab. 2016;101(12):4719–29. <https://pubmed.ncbi.nlm.nih.gov/27715326/>
63. Corbin JA, Bhaskar V, Goldfine ID, et al. **Improved glucose metabolism in vitro and in vivo by an allosteric monoclonal antibody that increases insulin receptor binding affinity.** PLoS One. 2014;9(2):e88684. <https://pubmed.ncbi.nlm.nih.gov/24533135/>
64. Johnson MB, De Franco E, et al. / **RZ358 (ersodetug) — insulin receptor antagonist antibody in congenital hyperinsulinism** (the development programme and early clinical work). Search: <https://pubmed.ncbi.nlm.nih.gov/?term=RZ358+OR+ersodetug+hyperinsulinism>
65. Müller D, Zimmering M, Roehr CC. **Should nifedipine be used to counter low blood sugar levels in children with persistent hyperinsulinaemic hypoglycaemia?** Arch Dis Child. 2004;89(1):83–5. <https://pubmed.ncbi.nlm.nih.gov/14709521/>

> Model dependency: entry 57 is the comparator for the reduction in glucose requirement on continuous glucagon infusion (35 % in the model).
> Entries 59 and 60 are the exendin(9-39) node and its PK. The **contrast** between entries 61 and 62 is the basis for writing sirolimus as
> "mechanism unknown, effect modest", and the model gives 19 % (stating that uncertainty at T4).
> Entries 63 and 64 ground the point that ersodetug works at the **receptor level**, that is independently of g (with the quantitative limits recorded at T3).
> Entry 65 is the clinical confirmation of the nifedipine failure, and the model explains that failure
> on **pharmacokinetic** grounds — "the mechanism is right, but the β-cell EC₅₀ is about 8 times the plasma Cmax" (P9).

---

## 9. Surgery — the intervention that reduces B and leaves g

66. Adzick NS, De Leon DD, States LJ, et al. **Surgical treatment of congenital hyperinsulinism: results from 500 pancreatectomies in neonates and children.** J Pediatr Surg. 2019;54(1):27–32. <https://pubmed.ncbi.nlm.nih.gov/30343978/>
67. Lord K, Dzata E, Snider KE, et al. **Clinical presentation and management of children with diffuse and focal hyperinsulinism: a review of 223 cases.** J Clin Endocrinol Metab. 2013;98(11):E1786–9. <https://pubmed.ncbi.nlm.nih.gov/24014812/>
68. Beltrand J, Caquard M, Arnoux JB, et al. **Glucose metabolism in 105 children and adolescents after pancreatectomy for congenital hyperinsulinism.** Diabetes Care. 2012;35(2):198–203. <https://pubmed.ncbi.nlm.nih.gov/22190679/>
69. Arya VB, Senniappan S, Demirbilek H, et al. **Pancreatic endocrine and exocrine function in children following near-total pancreatectomy for diffuse congenital hyperinsulinism.** PLoS One. 2014;9(5):e98054. <https://pubmed.ncbi.nlm.nih.gov/24840042/>
70. Lord K, Radcliffe J, Gallagher PR, et al. **High risk of diabetes and neurobehavioral deficits in individuals with surgically treated hyperinsulinism.** J Clin Endocrinol Metab. 2015;100(11):4133–9. <https://pubmed.ncbi.nlm.nih.gov/26327482/>

> Model dependency: **entries 68 and 70 are the comparators for P7.** In the model, surgery changes only `BMASS0` and
> does not touch `g`, so a single equation produces (i) persistent hypoglycaemia at a residual mass of 0.5,
> (ii) normal glucose at 0.2–0.3, and (iii) immediate diabetes at ≤0.1.
> The argument that growth reads the same residual mass again (2 % of 3.5 kg is 0.35 % at 20 kg)
> explains the "about half are diabetic by adolescence" of entries 68 and 70 **without any new assumption**.

---

## 10. Brain fuel — the reserve fuel that insulin deletes

71. Cremer JE. **Substrate utilization and brain development.** J Cereb Blood Flow Metab. 1982;2(4):394–407. <https://pubmed.ncbi.nlm.nih.gov/6754748/>
72. Nehlig A. **Brain uptake and metabolism of ketone bodies in animal models.** Prostaglandins Leukot Essent Fatty Acids. 2004;70(3):265–75. <https://pubmed.ncbi.nlm.nih.gov/14769485/>
73. Vannucci SJ, Simpson IA. **Developmental switch in brain nutrient transporter expression in the rat.** Am J Physiol Endocrinol Metab. 2003;285(5):E1127–34. <https://pubmed.ncbi.nlm.nih.gov/14534077/>
74. Pellerin L, Pellegri G, Bittar PG, et al. **Evidence supporting the existence of an activity-dependent astrocyte-neuron lactate shuttle.** Dev Neurosci. 1998;20(4-5):291–9. <https://pubmed.ncbi.nlm.nih.gov/9778565/>
75. Bougneres PF, Lemmel C, Ferré P, Bier DM. **Ketone body transport in the human neonate and infant.** J Clin Invest. 1986;77(1):42–8. <https://pubmed.ncbi.nlm.nih.gov/2867094/>
76. Kinnala A, Rikalainen H, Lapinleimu H, et al. **Cerebral magnetic resonance imaging and ultrasonography findings after neonatal hypoglycemia.** Pediatrics. 1999;103(4 Pt 1):724–9. <https://pubmed.ncbi.nlm.nih.gov/10103293/>
77. Burns CM, Rutherford MA, Boardman JP, Cowan FM. **Patterns of cerebral injury and neurodevelopmental outcomes after symptomatic neonatal hypoglycemia.** Pediatrics. 2008;122(1):65–74. <https://pubmed.ncbi.nlm.nih.gov/18595988/>
78. Menni F, de Lonlay P, Sevin C, et al. **Neurologic outcomes of 90 neonates and infants with persistent hyperinsulinemic hypoglycemia.** Pediatrics. 2001;107(3):476–9. <https://pubmed.ncbi.nlm.nih.gov/11230585/>
79. Avatapalle HB, Banerjee I, Shah S, et al. **Abnormal neurodevelopmental outcomes are common in children with transient congenital hyperinsulinism.** Front Endocrinol. 2013;4:60. <https://pubmed.ncbi.nlm.nih.gov/23730298/>
80. Ludwig A, Enke S, Heindorf J, et al. **Formal neurocognitive testing in 60 patients with congenital hyperinsulinism.** Horm Res Paediatr. 2018;89(1):1–6. <https://pubmed.ncbi.nlm.nih.gov/29131017/>

> Model dependency: **the whole of P5.** Entries 71, 72, and 75 are the basis for the neonatal brain actually oxidising ketones in bulk,
> and the point that the Km of the brain MCT (≈6 mM, entries 72 and 74) is the *rate-limiting step* is decisive.
> With Km set at 6 mM, ketones cover only 19 % of requirement at a BOHB of 2 mM, and as a result the
> permissible glucose in ketotic hypoglycaemia comes out at 25 mg/dL, close to the clinical observation (the 30–40 mg/dL range).
> Setting Km wrongly at 1.5 mM gave a permissible glucose of 16 mg/dL, at odds with the clinic, and
> that discrepancy forced the choice of Km — a **structural constraint**, not a calibration.
> Entries 76–80 are the outcome (NDD) nodes.

---

## 11. The dose-response separation of insulin (the IC₅₀ ordering) and normal neonatal metabolism

81. Rizza RA, Mandarino LJ, Gerich JE. **Dose-response characteristics for effects of insulin on production and utilization of glucose in man.** Am J Physiol. 1981;240(6):E630–9. <https://pubmed.ncbi.nlm.nih.gov/7018254/>
82. Nurjhan N, Campbell PJ, Kennedy FP, et al. **Insulin dose-response characteristics for suppression of glycerol release and conversion to glucose in humans.** Diabetes. 1986;35(12):1326–31. <https://pubmed.ncbi.nlm.nih.gov/3536590/>
83. Keller U, Schnell H, Sonnenberg GE, et al. **Role of insulin in the regulation of ketone body metabolism in man.** Diabetes Care. 1980;3(1):137–41. <https://pubmed.ncbi.nlm.nih.gov/6997936/>
84. Bier DM, Leake RD, Haymond MW, et al. **Measurement of "true" glucose production rates in infancy and childhood with 6,6-dideuteroglucose.** Diabetes. 1977;26(11):1016–23. <https://pubmed.ncbi.nlm.nih.gov/198425/>
85. Kalhan SC, Savin SM, Adam PA. **Measurement of glucose turnover in human newborn with glucose-1-¹³C.** J Clin Endocrinol Metab. 1976;43(3):704–7. <https://pubmed.ncbi.nlm.nih.gov/956353/>
86. Shelley HJ, Neligan GA. **Neonatal hypoglycaemia.** Br Med Bull. 1966;22(1):34–9. <https://pubmed.ncbi.nlm.nih.gov/5321967/>
87. Hume R, Burchell A, Williams FL, Koh DK. **Glucose homeostasis in the newborn.** Early Hum Dev. 2005;81(1):95–101. <https://pubmed.ncbi.nlm.nih.gov/15707720/>
88. Stanley CA, Rozance PJ, Thornton PS, et al. **Re-evaluating "transitional neonatal hypoglycemia": mechanism and implications for management.** J Pediatr. 2015;166(6):1520–5. <https://pubmed.ncbi.nlm.nih.gov/25819173/>

> Model dependency: **the basis for N8 (the IC₅₀ ordering) and for N6 and N7.** Entries 81–83 show that lipolysis and ketogenesis
> are suppressed at **far lower** insulin concentrations than hepatic glucose production, and
> that ordering automatically produces the conclusion that "CHI hypoglycaemia is always non-ketotic".
> Entries 84 and 85 give a neonatal glucose turnover of 4–6 mg/kg/min, entries 86 and 87 the hepatic glycogen store at term, and
> entry 88 offers the view that the transitional hypoglycaemia of the normal newborn is **the same shift in the K_ATP set point**.

---

## 12. Diagnostic tests — which term of the model each test reads

89. Finegold DN, Stanley CA, Baker L. **Glycemic response to glucagon during fasting hypoglycemia: an aid in the diagnosis of hyperinsulinism.** J Pediatr. 1980;96(2):257–9. <https://pubmed.ncbi.nlm.nih.gov/6766008/>
90. Ferrara C, Patel P, Becker S, et al. **Biomarkers of insulin for the diagnosis of hyperinsulinemic hypoglycemia in infants and children.** J Pediatr. 2016;168:212–9. <https://pubmed.ncbi.nlm.nih.gov/26490124/>
91. Palladino AA, Bennett MJ, Stanley CA. **Hyperinsulinism in infancy and childhood: when an insulin level is not always enough.** Clin Chem. 2008;54(2):256–63. <https://pubmed.ncbi.nlm.nih.gov/18156285/>
92. Hussain K, Hindmarsh P, Aynsley-Green A. **Neonates with symptomatic hyperinsulinemic hypoglycemia generate inappropriately low serum cortisol counterregulatory hormonal responses.** J Clin Endocrinol Metab. 2003;88(9):4342–7. <https://pubmed.ncbi.nlm.nih.gov/12970308/>

> Model dependency: **entry 89 is the direct counterpart of P6.** The glucagon stimulation test is diagnostic because
> insulin *keeps the glycogen from being used*, and the choice to make glycogen a
> **state variable** in the model automatically produces the discriminating power of this test (CHI +62 vs normal +16 mg/dL after a 20 h fast).
> Entry 91 is the insulin:C-peptide ratio node and entry 92 the counter-regulation failure node — in the model these are
> not separate lesions but a **consequence** of the α cell being suppressed by local insulin.

---

## 13. Hypothalamic K_ATP · hypoglycaemia awareness (the second organ of the same channel)

93. Miki T, Liss B, Minami K, et al. **ATP-sensitive K⁺ channels in the hypothalamus are essential for the maintenance of glucose homeostasis.** Nat Neurosci. 2001;4(5):507–12. <https://pubmed.ncbi.nlm.nih.gov/11319559/>
94. Evans ML, McCrimmon RJ, Flanagan DE, et al. **Hypothalamic ATP-sensitive K⁺ channels play a key role in sensing hypoglycemia and triggering counterregulatory epinephrine and glucagon responses.** Diabetes. 2004;53(10):2542–51. <https://pubmed.ncbi.nlm.nih.gov/15448088/>

> Model dependency: the dashed arrow `g → VMH` on the map. Because the same channel is also in the hypothalamus,
> in K_ATP-CHI the patient's **very organ for sensing hypoglycaemia is itself mutated**.
> The model does not quantify this term and leaves it as structure only (noted here for honesty).

---

## 14. Tools · methodology (QSP / mrgsolve)

95. Baron KT, Elmokadem A, et al. **mrgsolve: Simulate from ODE-Based Models.** <https://mrgsolve.org/>
96. Bergman RN, Ider YZ, Bowden CR, Cobelli C. **Quantitative estimation of insulin sensitivity.** Am J Physiol. 1979;236(6):E667–77. <https://pubmed.ncbi.nlm.nih.gov/443421/>
97. Dalla Man C, Rizza RA, Cobelli C. **Meal simulation model of the glucose-insulin system.** IEEE Trans Biomed Eng. 2007;54(10):1740–9. <https://pubmed.ncbi.nlm.nih.gov/17926672/>
98. Topp B, Promislow K, deVries G, et al. **A model of beta-cell mass, insulin, and glucose kinetics: pathways to diabetes.** J Theor Biol. 2000;206(4):605–19. <https://pubmed.ncbi.nlm.nih.gov/?term=Topp+model+beta-cell+mass+insulin+glucose+kinetics>
99. Fridlyand LE, Philipson LH. **Glucose sensing in the pancreatic beta cell: a computational systems analysis.** Theor Biol Med Model. 2010;7:15. <https://pubmed.ncbi.nlm.nih.gov/20497556/>
100. Chay TR, Keizer J. **Minimal model for membrane oscillations in the pancreatic beta-cell.** Biophys J. 1983;42(2):181–90. <https://pubmed.ncbi.nlm.nih.gov/6305437/>

> Model dependency: entries 96 and 97 are the lineage of the whole-body glucose-insulin compartment structure,
> entry 98 the lineage of the choice to make β-cell mass (B) a state variable, and
> entries 99 and 100 the lineage of the voltage divider and the `Popen → V_m → Ca` chain.
> What is new in this model is joining those two lineages through **a single parameter g** and
> classifying the drugs by their dependence on g (multiplicative / additive / independent).

---

## Appendix — Entries whose links are unstable

Entries 18, 48, 64, and 98 have been left with a search URL or an incomplete citation, either because the PMID is not certain
or because the material is at the development stage. Those references support the **structure** of the model (the trafficking defect, neonatal diazoxide PK,
the insulin receptor antibody, β-cell mass modelling), but not one **calibration number**
was taken from them.
