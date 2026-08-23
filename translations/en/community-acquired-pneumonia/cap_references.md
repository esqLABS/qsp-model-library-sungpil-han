# Community-Acquired Pneumonia (CAP) — References

The literature used for the structure, parameters and validation anchors of the QSP model
(`cap_qsp_model.dot`, `cap_mrgsolve_model.R`, `cap_shiny_app.R`). Every PMID was confirmed by
checking its title against PubMed E-utilities.

The **▶ Model link** item at the end of each section states which equation or parameter of the model
that literature supports.

---

## 1. Epidemiology · burden of disease · overview

1. Torres A, Cilloniz C, Niederman MS, et al. **Pneumonia.** *Nat Rev Dis Primers.* 2021;7:25.
   <https://pubmed.ncbi.nlm.nih.gov/33833230/>
2. Musher DM, Thorner AR. **Community-acquired pneumonia.** *N Engl J Med.* 2014;371:1619-1628.
   <https://pubmed.ncbi.nlm.nih.gov/25337751/>
3. Wunderink RG, Waterer GW. **Clinical practice. Community-acquired pneumonia.**
   *N Engl J Med.* 2014;370:543-551. <https://pubmed.ncbi.nlm.nih.gov/24499212/>
4. Prina E, Ranzani OT, Torres A. **Community-acquired pneumonia.** *Lancet.* 2015;386:1097-1108.
   <https://pubmed.ncbi.nlm.nih.gov/26277247/>
5. Jain S, Self WH, Wunderink RG, et al. **Community-acquired pneumonia requiring
   hospitalization among U.S. adults (EPIC study).** *N Engl J Med.* 2015;373:415-427.
   <https://pubmed.ncbi.nlm.nih.gov/26172429/>
6. Welte T, Torres A, Nathwani D. **Clinical and economic burden of community-acquired
   pneumonia among adults in Europe.** *Thorax.* 2012;67:71-79.
   <https://pubmed.ncbi.nlm.nih.gov/20729232/>
7. Ramirez JA, Wiemken TL, Peyrani P, et al. **Adults hospitalized with pneumonia in the
   United States: incidence, epidemiology, and mortality.** *Clin Infect Dis.* 2017;65:1806-1812.
   <https://pubmed.ncbi.nlm.nih.gov/29020164/>
8. Mizgerd JP. **Lung infection — a public health priority.** *PLoS Med.* 2006;3:e76.
   <https://pubmed.ncbi.nlm.nih.gov/16401173/>

**▶ Model link**: the range of the distributions of the initial inoculum `B0` and severity `SEV0`, and
the setting of the `HAZ0`/`HAZ_AGE` baseline hazard scale.

---

## 2. Clinical guidelines & severity scoring

9. Metlay JP, Waterer GW, Long AC, et al. **Diagnosis and treatment of adults with
   community-acquired pneumonia. An official ATS/IDSA clinical practice guideline.**
   *Am J Respir Crit Care Med.* 2019;200:e45-e67.
   <https://pubmed.ncbi.nlm.nih.gov/31573350/>
10. Mandell LA, Wunderink RG, Anzueto A, et al. **IDSA/ATS consensus guidelines on the
    management of community-acquired pneumonia in adults.** *Clin Infect Dis.* 2007;44 Suppl 2:S27-72.
    <https://pubmed.ncbi.nlm.nih.gov/17278083/>
11. Lim WS, van der Eerden MM, Laing R, et al. **Defining community acquired pneumonia
    severity on presentation to hospital (CURB-65).** *Thorax.* 2003;58:377-382.
    <https://pubmed.ncbi.nlm.nih.gov/12728155/>
12. Fine MJ, Auble TE, Yealy DM, et al. **A prediction rule to identify low-risk patients
    with community-acquired pneumonia (PSI).** *N Engl J Med.* 1997;336:243-250.
    <https://pubmed.ncbi.nlm.nih.gov/8995086/>
13. Singer M, Deutschman CS, Seymour CW, et al. **The Third International Consensus
    Definitions for Sepsis and Septic Shock (Sepsis-3).** *JAMA.* 2016;315:801-810.
    <https://pubmed.ncbi.nlm.nih.gov/26903338/>
14. Evans L, Rhodes A, Alhazzani W, et al. **Surviving Sepsis Campaign: international
    guidelines for management of sepsis and septic shock 2021.** *Crit Care Med.* 2021;49:e1063-e1143.
    <https://pubmed.ncbi.nlm.nih.gov/34643578/>
15. Salluh JIF, Souza-Dantas VC, Martin-Loeches I, et al. **Challenges for a broad
    international implementation of the current severe community-acquired pneumonia
    guidelines.** *Intensive Care Med.* 2024;50:526-538.
    <https://pubmed.ncbi.nlm.nih.gov/38546855/>

**▶ Model link**: implements as written the CURB-65 components in `[TABLE]` (confusion · `AKID` as a surrogate for urea · respiratory rate · blood pressure · age),
the SOFA sub-score cut-offs, and the Halm clinical stability criteria (`Clinically_stable`).

---

## 3. Pathogen biology & virulence factors

16. Weiser JN, Ferreira DM, Paton JC. ***Streptococcus pneumoniae*: transmission,
    colonization and invasion.** *Nat Rev Microbiol.* 2018;16:355-367.
    <https://pubmed.ncbi.nlm.nih.gov/29599457/>
17. Kadioglu A, Weiser JN, Paton JC, Andrew PW. **The role of *Streptococcus pneumoniae*
    virulence factors in host respiratory colonization and disease.** *Nat Rev Microbiol.*
    2008;6:288-301. <https://pubmed.ncbi.nlm.nih.gov/18340341/>
18. Loughran AJ, Orihuela CJ, Tuomanen EI. ***Streptococcus pneumoniae*: invasion and
    inflammation.** *Microbiol Spectr.* 2019;7. <https://pubmed.ncbi.nlm.nih.gov/30873934/>
19. Marriott HM, Mitchell TJ, Dockrell DH. **Pneumolysin: a double-edged sword during the
    host-pathogen interaction.** *Curr Mol Med.* 2008;8:497-509.
    <https://pubmed.ncbi.nlm.nih.gov/18781957/>
20. Bewley MA, Naughton M, Preston J, et al. **Pneumolysin activates macrophage lysosomal
    membrane permeabilization and executes apoptosis by distinct mechanisms.** *mBio.* 2014;5:e01710-14.
    <https://pubmed.ncbi.nlm.nih.gov/25293758/>
21. Anderson R, Nel JG, Feldman C. **Multifaceted role of pneumolysin in the pathogenesis
    of myocardial injury in community-acquired pneumonia.** *Int J Mol Sci.* 2018;19:1147.
    <https://pubmed.ncbi.nlm.nih.gov/29641429/>
22. Mizgerd JP. **Pathogenesis of severe pneumonia: advances and knowledge gaps.**
    *Curr Opin Pulm Med.* 2017;23:193-197. <https://pubmed.ncbi.nlm.nih.gov/28221171/>
23. McCullers JA. **Insights into the interaction between influenza virus and pneumococcus.**
    *Clin Microbiol Rev.* 2006;19:571-582. <https://pubmed.ncbi.nlm.nih.gov/16847087/>

**▶ Model link**: saturation of phagocytosis by the capsule (`KM_PHAG`), epithelial injury and TLR4/NLRP3 stimulation by pneumolysin
(`PLYFRAC`), the influenza NanA synergy (`EVIR_ADH`, `EVIR_BAR`).

---

## 4. Lysis-coupled PAMP liberation — the central term of this model

24. Nau R, Eiffert H. **Modulation of release of proinflammatory bacterial compounds by
    antibacterials: potential impact on course and outcome of bacterial infections.**
    *Clin Microbiol Rev.* 2002;15:95-110. <https://pubmed.ncbi.nlm.nih.gov/11781269/>
25. Spreer A, Kerstan H, Böttcher T, et al. **Reduced release of pneumolysin by
    *Streptococcus pneumoniae* in vitro and in vivo after treatment with nonbacteriolytic
    antibiotics.** *Antimicrob Agents Chemother.* 2003;47:2649-2654.
    <https://pubmed.ncbi.nlm.nih.gov/12878534/>
26. Lepper PM, Held TK, Schneider EM, et al. **Clinical implications of antibiotic-induced
    endotoxin release in septic shock.** *Intensive Care Med.* 2002;28:824-833.
    <https://pubmed.ncbi.nlm.nih.gov/12122518/>

**▶ Model link**: the whole `lysis_release = (YLYS_BL*kbl + YLYS_FQ*klvx + YLYS_MAC*kazi)*BE`
term of `dxdt_PAMP`. The basis for putting the lytic yield of a β-lactam at 1.0, a quinolone at 0.35 and
a macrolide at 0.12, and the only mechanism that produces the inflammatory surge and transient
clinical deterioration 6-24 hours after the first dose. Setting `YLYS_BL = 0` makes it disappear, so the structure is falsifiable.

---

## 5. Macrolide immunomodulation — a route to benefit without killing

27. Kovaleva A, Remmelts HHF, Rijkers GT, et al. **Immunomodulatory effects of macrolides
    during community-acquired pneumonia: a literature review.** *J Antimicrob Chemother.*
    2012;67:530-540. <https://pubmed.ncbi.nlm.nih.gov/22190609/>
28. Martin-Loeches I, Lisboa T, Rodriguez A, et al. **Combination antibiotic therapy with
    macrolides improves survival in intubated patients with community-acquired pneumonia.**
    *Intensive Care Med.* 2010;36:612-620. <https://pubmed.ncbi.nlm.nih.gov/19953222/>
29. Sligl WI, Asadi L, Eurich DT, et al. **Macrolides and mortality in critically ill
    patients with community-acquired pneumonia: a systematic review and meta-analysis.**
    *Crit Care Med.* 2014;42:420-432. <https://pubmed.ncbi.nlm.nih.gov/24158175/>
30. Restrepo MI, Anzueto A. **Macrolide therapy of pneumonia: is it necessary, and how does
    it help?** *Curr Opin Infect Dis.* 2016;29:212-217.
    <https://pubmed.ncbi.nlm.nih.gov/26836375/>

**▶ Model link**: `MAC_ANTI` (suppression of the cytokine amplifier, `IMAX_MAC`), `MAC_PLY` (inhibition of pneumolysin,
`IMAX_PLY`), `MAC_MIG` (inhibition of neutrophil migration, `IMAX_MIGM`) — all three are independent of
`kazi` (the bactericidal route). So most of the benefit remains even with `MIC_AZI = 64` (ermB resistance).

---

## 6. Antibiotic pharmacokinetics and pharmacodynamics

31. Craig WA. **Pharmacokinetic/pharmacodynamic parameters: rationale for antibacterial
    dosing of mice and men.** *Clin Infect Dis.* 1998;26:1-10.
    <https://pubmed.ncbi.nlm.nih.gov/9455502/>
32. Drusano GL. **Antimicrobial pharmacodynamics: critical interactions of 'bug and drug'.**
    *Nat Rev Microbiol.* 2004;2:289-300. <https://pubmed.ncbi.nlm.nih.gov/15031728/>
33. Craig WA, Andes D. **Pharmacokinetics and pharmacodynamics of antibiotics in otitis
    media.** *Pediatr Infect Dis J.* 1996;15:255-259.
    <https://pubmed.ncbi.nlm.nih.gov/8852915/>
34. Preston SL, Drusano GL, Berman AL, et al. **Pharmacodynamics of levofloxacin: a new
    paradigm for early clinical trials.** *JAMA.* 1998;279:125-129.
    <https://pubmed.ncbi.nlm.nih.gov/9440662/>
35. Drwiega EN, Rodvold KA. **Penetration of antibacterial agents into pulmonary epithelial
    lining fluid: an update.** *Clin Pharmacokinet.* 2022;61:17-46.
    <https://pubmed.ncbi.nlm.nih.gov/34651282/>
36. Kiem S, Schentag JJ. **Interpretation of epithelial lining fluid concentrations of
    antibiotics.** *Infect Chemother.* 2014;46:219-225.
    <https://pubmed.ncbi.nlm.nih.gov/25566401/>
37. Nightingale CH. **Pharmacokinetics and pharmacodynamics of newer macrolides.**
    *Pediatr Infect Dis J.* 1997;16:438-443. <https://pubmed.ncbi.nlm.nih.gov/9109156/>
38. Rapp RP, Kuhn R. **Comparison of five beta-lactam antibiotics against common nosocomial
    pathogens using the time above MIC at different creatinine clearances.**
    *Pharmacotherapy.* 2000;20:184S-190S. <https://pubmed.ncbi.nlm.nih.gov/10809349/>

**▶ Model link**: the `EMAX_*`/`FEC50_*`/`H_*` Emax-Hill structure and the ELF penetration ratios (`PEN_CEF` 0.40,
`PEN_AMX` 0.30, `PEN_LVX` 1.20, `KP_AZI_ELF` 22), and the `TAM` (fT>MIC) and `AUCF` (fAUC/MIC)
accumulating compartments. **Important**: these two indices are not parameters but the result of
integration. The difference between the β-lactam and the quinolone comes from one thing, the Hill coefficient (1.5 vs 2.6).

---

## 7. Corticosteroids — an intervention with two arms

39. Dequin PF, Meziani F, Quenot JP, et al. **Hydrocortisone in severe community-acquired
    pneumonia (CAPE-COD).** *N Engl J Med.* 2023;388:1931-1941.
    <https://pubmed.ncbi.nlm.nih.gov/36942789/>
40. Torres A, Sibila O, Ferrer M, et al. **Effect of corticosteroids on treatment failure
    among hospitalized patients with severe community-acquired pneumonia and high
    inflammatory response.** *JAMA.* 2015;313:677-686.
    <https://pubmed.ncbi.nlm.nih.gov/25688779/>
41. Blum CA, Nigro N, Briel M, et al. **Adjunct prednisone therapy for patients with
    community-acquired pneumonia: a multicentre, double-blind, randomised, placebo-
    controlled trial.** *Lancet.* 2015;385:1511-1518.
    <https://pubmed.ncbi.nlm.nih.gov/25608756/>
42. Confalonieri M, Urbino R, Potena A, et al. **Hydrocortisone infusion for severe
    community-acquired pneumonia: a preliminary randomized study.** *Am J Respir Crit Care Med.*
    2005;171:242-248. <https://pubmed.ncbi.nlm.nih.gov/15557131/>
43. Meduri GU, Shih MC, Bridges L, et al. **Low-dose methylprednisolone treatment in
    critically ill patients with severe community-acquired pneumonia (ESCAPe).**
    *Intensive Care Med.* 2022;48:1009-1023. <https://pubmed.ncbi.nlm.nih.gov/35723686/>
44. Saleem N, Kulkarni A, Snow TAC, et al. **Effect of corticosteroids on mortality and
    clinical cure in community-acquired pneumonia: a systematic review, meta-analysis,
    and trial sequential analysis.** *Chest.* 2023;163:484-497.
    <https://pubmed.ncbi.nlm.nih.gov/36087797/>
45. Stern A, Skalsky K, Avni T, et al. **Corticosteroids for pneumonia.**
    *Cochrane Database Syst Rev.* 2017;12:CD007720.
    <https://pubmed.ncbi.nlm.nih.gov/29236286/>
46. Feldman C, Anderson R. **Corticosteroids in the adjunctive therapy of community-acquired
    pneumonia: an appraisal of recent meta-analyses of clinical trials.**
    *J Thorac Dis.* 2016;8:E162-171. <https://pubmed.ncbi.nlm.nih.gov/27076965/>

**▶ Model link**: the two directions part company inside a single effect compartment, `HCE` — the
protective arm `GC_ANTI` (NF-κB transrepression, `IMAX_GC` 0.72) and the harmful arm `GC_PHAG` (inhibition of opsonic phagocytosis,
`IMAX_GCP` 0.45). Scenarios 6 and 7 use an identical steroid dose and differ only in MIC,
and what this model claims is that the **sign** of the net effect **inverts** there.
Setting `IMAX_GCP = 0` collapses that claim.

---

## 8. Biomarkers and treatment duration

47. Schuetz P, Wirz Y, Sager R, et al. **Effect of procalcitonin-guided antibiotic treatment
    on mortality in acute respiratory infections: a patient level meta-analysis.**
    *Lancet Infect Dis.* 2018;18:95-107. <https://pubmed.ncbi.nlm.nih.gov/29037960/>
48. Schuetz P, Christ-Crain M, Thomann R, et al. **Effect of procalcitonin-based guidelines
    vs standard guidelines on antibiotic use in lower respiratory tract infections
    (ProHOSP).** *JAMA.* 2009;302:1059-1066. <https://pubmed.ncbi.nlm.nih.gov/19738090/>
49. Chalmers JD, Singanayagam A, Hill AT. **C-reactive protein is an independent predictor
    of severity in community-acquired pneumonia.** *Am J Med.* 2008;121:219-225.
    <https://pubmed.ncbi.nlm.nih.gov/18328306/>
50. Uranga A, España PP, Bilbao A, et al. **Duration of antibiotic treatment in
    community-acquired pneumonia: a multicenter randomized clinical trial.**
    *JAMA Intern Med.* 2016;176:1257-1265. <https://pubmed.ncbi.nlm.nih.gov/27455166/>
51. Dinh A, Ropers J, Duran C, et al. **Discontinuing β-lactam treatment after 3 days for
    patients with community-acquired pneumonia in non-critical care wards (PTC): a
    double-blind, randomised, placebo-controlled, non-inferiority trial.**
    *Lancet.* 2021;397:1195-1203. <https://pubmed.ncbi.nlm.nih.gov/33773631/>
52. el Moussaoui R, de Borgie CAJM, van den Broek P, et al. **Effectiveness of discontinuing
    antibiotic treatment after three days versus eight days in mild to moderate-severe
    community acquired pneumonia.** *BMJ.* 2006;332:1355.
    <https://pubmed.ncbi.nlm.nih.gov/16763247/>

**▶ Model link**: why CRP (t½ 19 h, driven by IL-6) and PCT (t½ 24 h, driven by bacterial PAMP plus
interferon-mediated suppression `IPCT_VIR`) are separated onto different drivers. The safety of a short
course is structured so that it is decided not by `BE` but by the persister subpopulation `BP` (`KPERS`, `FTOL_BL`).

---

## 9. Outcomes · timing of treatment · long-term sequelae

53. Halm EA, Fine MJ, Marrie TJ, et al. **Time to clinical stability in patients hospitalized
    with community-acquired pneumonia: implications for practice guidelines.**
    *JAMA.* 1998;279:1452-1457. <https://pubmed.ncbi.nlm.nih.gov/9600479/>
54. Kumar A, Roberts D, Wood KE, et al. **Duration of hypotension before initiation of
    effective antimicrobial therapy is the critical determinant of survival in human septic
    shock.** *Crit Care Med.* 2006;34:1589-1596.
    <https://pubmed.ncbi.nlm.nih.gov/16625125/>
55. Houck PM, Bratzler DW. **Administration of first hospital antibiotics for
    community-acquired pneumonia: does timeliness affect outcomes?**
    *Curr Opin Infect Dis.* 2005;18:151-156. <https://pubmed.ncbi.nlm.nih.gov/15735420/>
56. Corrales-Medina VF, Alvarez KN, Weissfeld LA, et al. **Association between
    hospitalization for pneumonia and subsequent risk of cardiovascular disease.**
    *JAMA.* 2015;313:264-274. <https://pubmed.ncbi.nlm.nih.gov/25602997/>
57. Corrales-Medina VF, Taljaard M, Fine MJ, et al. **Risk stratification for cardiac
    complications in patients hospitalized for community-acquired pneumonia.**
    *Mayo Clin Proc.* 2014;89:60-68. <https://pubmed.ncbi.nlm.nih.gov/24388023/>
58. Rudd KE, Johnson SC, Agesa KM, et al. **Global, regional, and national sepsis incidence
    and mortality, 1990-2017: analysis for the Global Burden of Disease Study.**
    *Lancet.* 2020;395:200-211. <https://pubmed.ncbi.nlm.nih.gov/31954465/>
59. Piqueras M, et al. **Persistent systemic interleukin-6 elevation after community-acquired
    pneumonia is associated with one-year outcomes.** *Respir Res.* 2025.
    <https://pubmed.ncbi.nlm.nih.gov/41239435/>

**▶ Model link**: the `dxdt_CHZ` cumulative hazard function (weighted on SOFA · age · bacteraemic log CFU) and the
door-to-antibiotic delay geometry of scenario 9 — hazard accumulated during the delay is never recovered afterwards.

---

## 10. Lung injury · barrier · gas exchange · resolution

60. Matthay MA, Zemans RL, Zimmerman GA, et al. **Acute respiratory distress syndrome.**
    *Nat Rev Dis Primers.* 2019;5:18. <https://pubmed.ncbi.nlm.nih.gov/30872586/>
61. Ware LB, Matthay MA. **Alveolar fluid clearance is impaired in the majority of patients
    with acute lung injury and the acute respiratory distress syndrome.**
    *Am J Respir Crit Care Med.* 2001;163:1376-1383.
    <https://pubmed.ncbi.nlm.nih.gov/11371404/>
62. Roux J, McNicholas CM, Carles M, et al. **IL-8 inhibits cAMP-stimulated alveolar
    epithelial fluid transport via a GRK2/PI3K-dependent mechanism.** *FASEB J.*
    2013;27:1095-1106. <https://pubmed.ncbi.nlm.nih.gov/23221335/>
63. Herold S, Mayer K, Lohmeyer J. **Acute lung injury: how macrophages orchestrate
    resolution of inflammation and tissue repair.** *Front Immunol.* 2011;2:65.
    <https://pubmed.ncbi.nlm.nih.gov/22566854/>
64. Serhan CN. **Pro-resolving lipid mediators are leads for resolution physiology.**
    *Nature.* 2014;510:92-101. <https://pubmed.ncbi.nlm.nih.gov/24899309/>
65. Aberdein JD, Cole J, Bewley MA, et al. **Alveolar macrophages in pulmonary host defence
    — the unrecognised role of apoptosis as a mechanism of intracellular bacterial killing.**
    *Clin Exp Immunol.* 2013;174:193-202. <https://pubmed.ncbi.nlm.nih.gov/23841514/>
66. Dockrell DH, Whyte MKB, Mitchell TJ. **Pneumococcal pneumonia: mechanisms of infection
    and resolution.** *Chest.* 2012;142:482-491.
    <https://pubmed.ncbi.nlm.nih.gov/22871758/>
67. Grudzinska F, et al. **Hospitalised older adults with community-acquired pneumonia and
    sepsis have dysregulated neutrophil function.** *Thorax.* 2025.
    <https://pubmed.ncbi.nlm.nih.gov/39689942/>
68. Mohanty T, Fisher J, Bakochi A, et al. **Neutrophil extracellular traps in the central
    nervous system hinder bacterial clearance during pneumococcal meningitis.**
    *Nat Commun.* 2019;10:1667. <https://pubmed.ncbi.nlm.nih.gov/30971685/>
69. Archer SL, Dunham-Snary KJ, et al. **Hypoxic pulmonary vasoconstriction: an important
    component of the homeostatic oxygen sensing system.** *Physiol Res.* 2024;73:S493-S510.
    <https://pubmed.ncbi.nlm.nih.gov/39589299/>

**▶ Model link**: `IAFC_IL8` (inhibition of alveolar fluid clearance by IL-8, ref 62), `KS_SPM`/`ESPM`
(efferocytosis → SPM → damping of the amplifier, ref 63-64), the basis for `KPHAG_AM` being higher than
for the neutrophil (ref 65), and `HPV_MAX`/`IHPV_CYT` (hypoxic pulmonary vasoconstriction and its inflammatory blunting, ref 69).

---

## 11. Prevention

70. Bonten MJM, Huijts SM, Bolkenbaas M, et al. **Polysaccharide conjugate vaccine against
    pneumococcal pneumonia in adults (CAPiTA).** *N Engl J Med.* 2015;372:1114-1125.
    <https://pubmed.ncbi.nlm.nih.gov/25785969/>
71. Kobayashi M, Leidner AJ, Gierke R, et al. **Use of 21-valent pneumococcal conjugate
    vaccine among U.S. adults: recommendations of the Advisory Committee on Immunization
    Practices.** *MMWR Morb Mortal Wkly Rep.* 2024;73:793-798.
    <https://pubmed.ncbi.nlm.nih.gov/39264843/>

**▶ Model link**: the `VACC` covariate → a 1.45-fold rise in opsonic activity, `OPS_VACC`.

---

## Appendix: summary table of parameter provenance

| Parameter group | Value | Source |
|---|---|---|
| Ceftriaxone CL 1.0 L/h, V1 8 L, fu 0.10, ELF 0.40 | 2 g q24h → Cmax ~150, trough ~10 mg/L | 35, 36, 38 |
| Amoxicillin ka 1.5/h, F 0.80, CL 20 L/h, fu 0.82 | 1 g q8h | 33, 35 |
| Azithromycin F 0.37, Vc 450 L, deep compartment 2200 L, ELF factor 22 | serum ~0.4 → ELF 1-3 mg/L | 37, 35 |
| Levofloxacin CL 8 L/h, V 90 L, fu 0.70, ELF 1.2 | 750 mg → AUC24 ~100, fAUC/MIC ~70 | 34, 35 |
| Hydrocortisone CL 18 L/h, V 35 L, ke0 0.2/h | continuous infusion of 200 mg/day | 39, 42 |
| CRP t½ 19 h · PCT t½ 24 h | biomarker decay | 47-49 |
| Lytic yield β-lactam 1.0 / FQ 0.35 / macrolide 0.12 | PAMP release by antibiotic class | 24-26 |
| Steroid protective arm 0.72 / harmful arm 0.45 | reproducing the direction of CAPE-COD | 39-46 |
| Median ~3 days to clinical stability (on effective treatment) | validation target | 53 |
| The "crisis" at 7-9 days of the untreated natural history | validation target | 2, 3 |

---

## Disclaimer

This model and this literature compilation are **for education and research**. The parameters are
approximations based on the published literature and have not been fitted to patient-level data. They cannot be used for clinical decision-making, prescribing or regulatory submission.
