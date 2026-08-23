# Prosthetic Joint Infection (PJI) — References

> Implant-associated *Staphylococcus aureus* osteomyelitis — the evidence base for the QSP model
>
> **How the PMIDs were confirmed**: the 68 PMIDs below were all looked up with NCBI E-utilities (`esearch` + `esummary`)
> and checked against the title, journal, year, and first author. They are not numbers written from memory, and
> the journal and year of each entry are exactly as PubMed returned them.

---

## 0. The four pillars and where they come from

| Pillar | Quantitative claim | Primary evidence |
|------|-------------|----------|
| ① A foreign body shifts the ID50 | with a foreign body present, the inoculum needed to establish infection falls by 10⁴–10⁵ fold | Elek & Conen 1957 [13499821] · Zimmerli 1982 [7119479] |
| ② Killing is C_bone,free / MBEC | the free bone concentration achievable by systemic dosing is 0.004–0.28 of the MBEC | Ceri 1999 [10325322] · Landersdorfer 2009 [19271782] · Thabit 2019 [30772469] |
| ③ Surgery is manipulation of the mutant supply | P(a pre-existing rpoB mutant) = 1 − e^(−μN), μ ≈ 10⁻⁸ | O'Neill 2001 [11328777] · Achermann 2013 [22987291] · Drlica 2007 [17278059] |
| ④ Biofilm maturation is a clock | the success rate of DAIR is a function of the duration of symptoms | Byren 2009 [19336454] · Lora-Tamayo 2013 [22942204] · Nishitani 2015 [25820925] |

---

## 1. Reviews and disease overview

1. Zimmerli W, Trampuz A, Ochsner PE. **Prosthetic-joint infections.** *N Engl J Med.* 2004. — the standard review of the field. The prototype for the pathogenesis, diagnosis, and treatment algorithm.
   <https://pubmed.ncbi.nlm.nih.gov/15483283/>
2. Tande AJ, Patel R. **Prosthetic joint infection.** *Clin Microbiol Rev.* 2014.
   <https://pubmed.ncbi.nlm.nih.gov/24696437/>
3. Kavanagh N, et al. **Staphylococcal Osteomyelitis: Disease Progression, Treatment Challenges, and Future Directions.** *Clin Microbiol Rev.* 2018.
   <https://pubmed.ncbi.nlm.nih.gov/29444953/>
4. Masters EA, et al. **Skeletal infections: microbial pathogenesis, immunity and clinical management.** *Nat Rev Microbiol.* 2022.
   <https://pubmed.ncbi.nlm.nih.gov/35169289/>
5. Masters EA, et al. **Evolving concepts in bone infection: redefining "biofilm", "acute vs. chronic osteomyelitis", "the immune proteome" and "local antibiotic therapy".** *Bone Res.* 2019.
   <https://pubmed.ncbi.nlm.nih.gov/31646012/>
6. Lew DP, Waldvogel FA. **Osteomyelitis.** *Lancet.* 2004.
   <https://pubmed.ncbi.nlm.nih.gov/15276398/>
7. Izakovicova P, Borens O, Trampuz A. **Periprosthetic joint infection: current concepts and outlook.** *EFORT Open Rev.* 2019.
   <https://pubmed.ncbi.nlm.nih.gov/31423332/>
8. Zimmerli W, Sendi P. **Orthopaedic biofilm infections.** *APMIS.* 2017.
   <https://pubmed.ncbi.nlm.nih.gov/28407423/>

---

## 2. Pillar ① — the foreign body and the infective dose

9. Elek SD, Conen PE. **The virulence of *Staphylococcus pyogenes* for man; a study of the problems of wound infection.** *Br J Exp Pathol.* 1957. — the classic experiment showing a 10⁴–10⁵ fold difference in the inoculum that produces a lesion in human skin with and without a silk suture present. What the model's combination of `IMPL`/`FBFACT`/`KATT` sets out to reproduce.
   <https://pubmed.ncbi.nlm.nih.gov/13499821/>
10. Zimmerli W, Waldvogel FA, Vaudaux P, Nydegger UE. **Pathogenesis of foreign body infection: description and characteristics of an animal model.** *J Infect Dis.* 1982. — the tissue cage model. The experimental basis for the local neutrophil functional defect (the model's `FBFACT`).
    <https://pubmed.ncbi.nlm.nih.gov/7119479/>
11. Thurlow LR, et al. **Staphylococcus aureus biofilms prevent macrophage phagocytosis and attenuate inflammation in vivo.** *J Immunol.* 2011. — frustrated phagocytosis (`FRUST`).
    <https://pubmed.ncbi.nlm.nih.gov/21525381/>
12. Heim CE, et al. **Myeloid-derived suppressor cells contribute to Staphylococcus aureus orthopedic biofilm infection.** *J Immunol.* 2014. — the model's `MDSC` and `EMDSC` terms.
    <https://pubmed.ncbi.nlm.nih.gov/24646737/>
13. Heim CE, et al. **Interleukin-10 production by myeloid-derived suppressor cells contributes to bacterial persistence during Staphylococcus aureus orthopedic biofilm infection.** *J Leukoc Biol.* 2015.
    <https://pubmed.ncbi.nlm.nih.gov/26232453/>

---

## 3. Biofilm biology

14. Costerton JW, Stewart PS, Greenberg EP. **Bacterial biofilms: a common cause of persistent infections.** *Science.* 1999.
    <https://pubmed.ncbi.nlm.nih.gov/10334980/>
15. Stewart PS, Costerton JW. **Antibiotic resistance of bacteria in biofilms.** *Lancet.* 2001. — the triple mechanism of retarded diffusion, metabolic arrest, and persisters.
    <https://pubmed.ncbi.nlm.nih.gov/11463434/>
16. Hall CW, Mah TF. **Molecular mechanisms of biofilm-based antibiotic resistance and tolerance in pathogenic bacteria.** *FEMS Microbiol Rev.* 2017.
    <https://pubmed.ncbi.nlm.nih.gov/28369412/>
17. Moormeier DE, Bayles KW. **Staphylococcus aureus biofilm: a complex developmental organism.** *Mol Microbiol.* 2017. — the developmental stages of attachment → proliferation → dispersal (the model's `KATT`/`EPS`/`KDISP`).
    <https://pubmed.ncbi.nlm.nih.gov/28142193/>
18. McConoughey SJ, et al. **Biofilms in periprosthetic orthopedic infections.** *Future Microbiol.* 2014.
    <https://pubmed.ncbi.nlm.nih.gov/25302955/>
19. Nishitani K, et al. **Quantifying the natural history of biofilm formation in vivo during the establishment of chronic implant-associated Staphylococcus aureus osteomyelitis in mice.** *J Orthop Res.* 2015. — the experimental counterpart of the model's EPS maturation time constant (~10 days).
    <https://pubmed.ncbi.nlm.nih.gov/25820925/>
20. Foster TJ, Geoghegan JA, Ganesh VK, Höök M. **Adhesion, invasion and evasion: the many functions of the surface proteins of Staphylococcus aureus.** *Nat Rev Microbiol.* 2014. — MSCRAMM(FnBPA/ClfA/SdrC-E).
    <https://pubmed.ncbi.nlm.nih.gov/24336184/>
21. Geoghegan JA, Foster TJ. **Cell Wall-Anchored Surface Proteins of Staphylococcus aureus: Many Proteins, Multiple Functions.** *Curr Top Microbiol Immunol.* 2017.
    <https://pubmed.ncbi.nlm.nih.gov/26667044/>
22. Otto M. **Phenol-soluble modulins.** *Int J Med Microbiol.* 2014. — the PSMs: matrix amyloid, neutrophil lysis, dispersal.
    <https://pubmed.ncbi.nlm.nih.gov/24447915/>
23. Cheung GY, Joo HS, Chatterjee SS, Otto M. **Phenol-soluble modulins — critical determinants of staphylococcal virulence.** *FEMS Microbiol Rev.* 2014.
    <https://pubmed.ncbi.nlm.nih.gov/24372362/>

---

## 4. Phenotypic subpopulations — persisters, SCVs, the intracellular reservoir

24. Lewis K. **Persister cells.** *Annu Rev Microbiol.* 2010. — the model's `NPER`, `KPER`, and `KWAKE`.
    <https://pubmed.ncbi.nlm.nih.gov/20528688/>
25. Proctor RA, et al. **Small colony variants: a pathogenic form of bacteria that facilitates persistent and recurrent infections.** *Nat Rev Microbiol.* 2006. — `NSCV`.
    <https://pubmed.ncbi.nlm.nih.gov/16541137/>
26. de Mesy Bentley KL, et al. **Evidence of Staphylococcus aureus Deformation, Proliferation, and Migration in Canaliculi of Live Cortical Bone in Murine Models of Osteomyelitis.** *J Bone Miner Res.* 2017. — the canalicular reservoir: the basis for why the model's `NIC` is not readily killed by an antibiotic.
    <https://pubmed.ncbi.nlm.nih.gov/27933662/>
27. Ellington JK, et al. **Intracellular Staphylococcus aureus and antibiotic resistance: implications for treatment of staphylococcal osteomyelitis.** *J Orthop Res.* 2006. — the difference in activity between antibiotics against intracellular organisms (`CAR`, `PHI*`).
    <https://pubmed.ncbi.nlm.nih.gov/16419973/>
28. Josse J, Velard F, Gangloff SC. **Staphylococcus aureus vs. Osteoblast: Relationship and Consequences in Osteomyelitis.** *Front Cell Infect Microbiol.* 2015.
    <https://pubmed.ncbi.nlm.nih.gov/26636047/>

---

## 5. Pillar ② — MBEC, bone penetration, PK/PD (the C_bone,free / MBEC ratio)

29. Ceri H, et al. **The Calgary Biofilm Device: new technology for rapid determination of antibiotic susceptibilities of bacterial biofilms.** *J Clin Microbiol.* 1999. — the source of the MBEC concept and of how it is measured.
    <https://pubmed.ncbi.nlm.nih.gov/10325322/>
30. Landersdorfer CB, Bulitta JB, Kinzig M, Holzgrabe U, Sörgel F. **Penetration of antibacterials into bone: pharmacokinetic, pharmacodynamic and bioanalytical considerations.** *Clin Pharmacokinet.* 2009. — the primary source of the model's `PEN*` values.
    <https://pubmed.ncbi.nlm.nih.gov/19271782/>
31. Thabit AK, et al. **Antibiotic penetration into bone and joints: An updated review.** *Int J Infect Dis.* 2019.
    <https://pubmed.ncbi.nlm.nih.gov/30772469/>
32. Drusano GL. **Antimicrobial pharmacodynamics: critical interactions of 'bug and drug'.** *Nat Rev Microbiol.* 2004.
    <https://pubmed.ncbi.nlm.nih.gov/15031728/>
33. Post V, et al. **Vancomycin displays time-dependent eradication of mature Staphylococcus aureus biofilms.** *J Orthop Res.* 2017. — the observation that vancomycin's eradication of biofilm depends on the duration of exposure rather than on concentration (the model's combination of a low `EMXVAN` with a high `MBCVAN`).
    <https://pubmed.ncbi.nlm.nih.gov/27175462/>
34. Widmer AF, Frei R, Rajacic Z, Zimmerli W. **Correlation between in vivo and in vitro efficacy of antimicrobial agents against foreign body infections.** *J Infect Dis.* 1990. — direct evidence that in foreign body infection the in vitro MIC does not predict in vivo efficacy.
    <https://pubmed.ncbi.nlm.nih.gov/2355207/>
35. Rybak MJ, et al. **Therapeutic monitoring of vancomycin for serious methicillin-resistant Staphylococcus aureus infections: A revised consensus guideline and review.** *Am J Health Syst Pharm.* 2020. — the AUC24 target of 400–600.
    <https://pubmed.ncbi.nlm.nih.gov/32191793/>

---

## 6. Pillar ③ — rifampicin and the mutant supply

36. Zimmerli W, Widmer AF, Blatter M, Frei R, Ochsner PE. **Role of rifampin for treatment of orthopedic implant-related staphylococcal infections: a randomized controlled trial.** *JAMA.* 1998. — ciprofloxacin + rifampicin 12/12 (100%) against ciprofloxacin alone 7/12 (58%). The calibration reference point for scenarios 05/06 of the model.
    <https://pubmed.ncbi.nlm.nih.gov/9605897/>
37. Widmer AF, Gaechter A, Ochsner PE, Zimmerli W. **Antimicrobial treatment of orthopedic implant-related infections with rifampin combinations.** *Clin Infect Dis.* 1992.
    <https://pubmed.ncbi.nlm.nih.gov/1623081/>
38. Zimmerli W, Sendi P. **Role of Rifampin against Staphylococcal Biofilm Infections In Vitro, in Animal Models, and in Orthopedic-Device-Related Infections.** *Antimicrob Agents Chemother.* 2019.
    <https://pubmed.ncbi.nlm.nih.gov/30455229/>
39. Sendi P, Zimmerli W. **The use of rifampin in staphylococcal orthopaedic-device-related infections.** *Clin Microbiol Infect.* 2017.
    <https://pubmed.ncbi.nlm.nih.gov/27746393/>
40. Sendi P, Zimmerli W. **Antimicrobial treatment concepts for orthopaedic device-related infection.** *Clin Microbiol Infect.* 2012.
    <https://pubmed.ncbi.nlm.nih.gov/23046277/>
41. O'Neill AJ, Cove JH, Chopra I. **Mutation frequencies for resistance to fusidic acid and rifampicin in Staphylococcus aureus.** *J Antimicrob Chemother.* 2001. — the source of μ_rpoB ≈ 10⁻⁷–10⁻⁸. The model's `MURPOB`.
    <https://pubmed.ncbi.nlm.nih.gov/11328777/>
42. Aubry-Damon H, Soussy CJ, Courvalin P. **Characterization of mutations in the rpoB gene that confer rifampin resistance in Staphylococcus aureus.** *Antimicrob Agents Chemother.* 1998.
    <https://pubmed.ncbi.nlm.nih.gov/9756760/>
43. Achermann Y, et al. **Factors associated with rifampin resistance in staphylococcal periprosthetic joint infections (PJI): a matched case-control study.** *Infection.* 2013. — the conditions under which rifampicin resistance actually emerges in the clinic (inadequate debridement, exposure as monotherapy).
    <https://pubmed.ncbi.nlm.nih.gov/22987291/>
44. Drlica K, Zhao X. **Mutant selection window hypothesis updated.** *Clin Infect Dis.* 2007. — the MSW concept (the model's `RES_MSW` node).
    <https://pubmed.ncbi.nlm.nih.gov/17278059/>

---

## 7. Pillar ④ — surgical strategy and the maturity clock

45. Byren I, et al. **One hundred and twelve infected arthroplasties treated with 'DAIR' (debridement, antibiotics and implant retention): antibiotic duration and outcome.** *J Antimicrob Chemother.* 2009.
    <https://pubmed.ncbi.nlm.nih.gov/19336454/>
46. Lora-Tamayo J, et al. **A large multicenter study of methicillin-susceptible and methicillin-resistant Staphylococcus aureus prosthetic joint infections managed with implant retention.** *Clin Infect Dis.* 2013.
    <https://pubmed.ncbi.nlm.nih.gov/22942204/>
47. Trebse R, Pisot V, Trampuz A. **Treatment of infected retained implants.** *J Bone Joint Surg Br.* 2005.
    <https://pubmed.ncbi.nlm.nih.gov/15736752/>
48. Bejon P, et al. **Two-stage revision for prosthetic joint infection: predictors of outcome and the role of reimplantation microbiology.** *J Antimicrob Chemother.* 2010.
    <https://pubmed.ncbi.nlm.nih.gov/20053693/>
49. Kunutsor SK, Whitehouse MR, Blom AW, Beswick AD. **Re-Infection Outcomes following One- and Two-Stage Surgical Revision of Infected Hip Prosthesis: A Systematic Review and Meta-Analysis.** *PLoS One.* 2015.
    <https://pubmed.ncbi.nlm.nih.gov/26407003/>
50. Bertazzoni Minelli E, Benini A, Magnan B, Bartolozzi P. **Release of gentamicin and vancomycin from temporary human hip spacers in two-stage revision of infected arthroplasty.** *J Antimicrob Chemother.* 2004. — the basis for the model's two-pool spacer elution parameters (`KFAST`/`KSLOW`/`SPCF0`/`SPCS0`).
    <https://pubmed.ncbi.nlm.nih.gov/14688051/>

---

## 8. Clinical trials of antimicrobial regimens

51. Li HK, et al. (OVIVA). **Oral versus Intravenous Antibiotics for Bone and Joint Infection.** *N Engl J Med.* 2019. — oral was non-inferior to intravenous (failure 13.2% against 14.6%). Scenarios 11/12 of the model: what decides is not the route but the bone exposure.
    <https://pubmed.ncbi.nlm.nih.gov/30699315/>
52. Bernard L, et al. (DATIPO). **Antibiotic Therapy for 6 or 12 Weeks for Prosthetic Joint Infection.** *N Engl J Med.* 2021. — 6 weeks was **not** non-inferior to 12 weeks. The result scenarios 13/14 of the model reproduce.
    <https://pubmed.ncbi.nlm.nih.gov/34042388/>
53. Byren I, et al. **Randomized controlled trial of the safety and efficacy of Daptomycin versus standard-of-care therapy for management of patients with osteomyelitis associated with prosthetic devices undergoing two-stage revision arthroplasty.** *Antimicrob Agents Chemother.* 2012.
    <https://pubmed.ncbi.nlm.nih.gov/22908174/>
54. John AK, et al. **Efficacy of daptomycin in implant-associated infection due to methicillin-resistant Staphylococcus aureus: importance of combination with rifampin.** *Antimicrob Agents Chemother.* 2009.
    <https://pubmed.ncbi.nlm.nih.gov/19364845/>
55. Saleh-Mghir A, et al. **Adjunctive rifampin is crucial to optimizing daptomycin efficacy against rabbit prosthetic joint infection due to methicillin-resistant Staphylococcus aureus.** *Antimicrob Agents Chemother.* 2011.
    <https://pubmed.ncbi.nlm.nih.gov/21825285/>
56. Baldoni D, Haschke M, Rajacic Z, Zimmerli W, Trampuz A. **Linezolid alone or combined with rifampin against methicillin-resistant Staphylococcus aureus in experimental foreign-body infection.** *Antimicrob Agents Chemother.* 2009.
    <https://pubmed.ncbi.nlm.nih.gov/19075065/>

---

## 9. Diagnostic criteria and biomarkers

57. Parvizi J, Tan TL, Goswami K, et al. **The 2018 Definition of Periprosthetic Hip and Knee Infection: An Evidence-Based and Validated Criteria.** *J Arthroplasty.* 2018. — the model's `ICMSC` formula.
    <https://pubmed.ncbi.nlm.nih.gov/29551303/>
58. Osmon DR, et al. **Diagnosis and management of prosthetic joint infection: clinical practice guidelines by the Infectious Diseases Society of America.** *Clin Infect Dis.* 2013.
    <https://pubmed.ncbi.nlm.nih.gov/23223583/>
59. Trampuz A, et al. **Sonication of removed hip and knee prostheses for diagnosis of infection.** *N Engl J Med.* 2007. — direct evidence that the biofilm has to be dislodged before it will culture (the model's `NB` → `DX_SONIC` pathway).
    <https://pubmed.ncbi.nlm.nih.gov/17699815/>
60. Deirmengian C, et al. **Diagnosing periprosthetic joint infection: has the era of the biomarker arrived?** *Clin Orthop Relat Res.* 2014. — α-defensin.
    <https://pubmed.ncbi.nlm.nih.gov/24590839/>
61. Renz N, et al. **Alpha Defensin Lateral Flow Test for Diagnosis of Periprosthetic Joint Infection: Not a Screening but a Confirmatory Test.** *J Bone Joint Surg Am.* 2018.
    <https://pubmed.ncbi.nlm.nih.gov/29715222/>

---

## 10. Bone destruction

62. Nair SP, Meghji S, Wilson M, et al. **Bacterially induced bone destruction: mechanisms and misconceptions.** *Infect Immun.* 1996.
    <https://pubmed.ncbi.nlm.nih.gov/8698454/>
63. Boyce BF, Xing L. **Functions of RANKL/RANK/OPG in bone modeling and remodeling.** *Arch Biochem Biophys.* 2008. — the model's `RANKL`/`OPG`/`OCL` block.
    <https://pubmed.ncbi.nlm.nih.gov/18395508/>

---

## 11. Drug toxicity and interactions

64. Lodise TP, et al. **Relationship between initial vancomycin concentration-time profile and nephrotoxicity among hospitalized patients.** *Clin Infect Dis.* 2009. — the model's `AUCTHR` and `ENEPH`.
    <https://pubmed.ncbi.nlm.nih.gov/19586413/>
65. Boak LM, et al. **Clinical population pharmacokinetics and toxicodynamics of linezolid.** *Antimicrob Agents Chemother.* 2014. — the exposure-response of linezolid thrombocytopenia (`EPLT`/`KPLT`).
    <https://pubmed.ncbi.nlm.nih.gov/24514086/>
66. Gandelman K, et al. **Unexpected effect of rifampin on the pharmacokinetics of linezolid: in silico and in vitro approaches to explain its mechanism.** *J Clin Pharmacol.* 2011. — a ~32% fall in linezolid AUC caused by rifampicin. The model's `RIFLZD` term.
    <https://pubmed.ncbi.nlm.nih.gov/20371736/>

---

## 12. Burden and cost

67. Kurtz SM, Lau E, Watson H, Schmier JK, Parvizi J. **Economic burden of periprosthetic joint infection in the United States.** *J Arthroplasty.* 2012.
    <https://pubmed.ncbi.nlm.nih.gov/22554729/>
68. Kurtz SM, et al. **Hospital Costs for Unsuccessful Two-Stage Revisions for Periprosthetic Joint Infection.** *J Arthroplasty.* 2022.
    <https://pubmed.ncbi.nlm.nih.gov/34763048/>

---

## Appendix: Where the model parameters came from in the literature

| Model parameter | Value | Evidence |
|---------------|-----|------|
| `PENVAN` / `PENRIF` / `PENLVX` / `PENDAP` / `PENLZD` / `PENCFZ` | 0.20 / 0.50 / 0.50 / 0.15 / 0.45 / 0.20 | [19271782] · [30772469] |
| `MBCVAN` … `MBCRIF` | 512 … 1.0 mg/L | [10325322] · [27175462] · [30455229] |
| `MURPOB` | 1×10⁻⁸ /division | [11328777] |
| `FBFACT` (the foreign body immune defect) | 0.12 | [7119479] · [13499821] |
| `FRUST` (frustrated phagocytosis) | 0.03 | [21525381] |
| `KSMDSC` / `EMDSC` | 0.03 /h, 0.70 | [24646737] · [26232453] |
| `KFAST` / `KSLOW` / `SPCF0` / `SPCS0` | 0.05 /h, 0.0015 /h, 1200 mg, 700 mg | [14688051] |
| `AUCTHR` / `ENEPH` | 600 mg·h/L, 0.55 | [19586413] · [32191793] |
| `EPLT` / `KPLT` | 0.60, 9 mg/L | [24514086] |
| `RIFLZD` | 0.50 (LZD AUC −32%) | [20371736] |
| Validation reference point (scenarios 05/06) | RIF combination 100% against monotherapy 58% | [9605897] |
| Validation reference point (scenarios 11/12) | oral non-inferior | [30699315] |
| Validation reference point (scenarios 13/14) | 6 weeks ≠ 12 weeks | [34042388] |

---

*This reference list is part of the documentation of a QSP model built for educational and research purposes. Do not use it directly in clinical decision-making.*
