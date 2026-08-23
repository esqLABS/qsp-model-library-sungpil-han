# Influenza A — References

Every entry in this list is **bibliographic information looked up and confirmed** through the PubMed
E-utilities, and the author, year, journal, volume and pages were taken verbatim from the PubMed
record. No PMID written from memory is included.

The heading of each section states **which part of the model** those papers support.
Provenance is marked on the `$PARAM` block of `flu_mrgsolve_model.R` as `[LIT]`
(a literature value) / `[CAL]` (calibrated in this model) / `[ASM]` (a structural assumption), and
`[CAL]` and `[ASM]` mean the value is not a literature value.

> **Only two sources were actually used in the calibration.** The three endpoints of CAPSTONE-1
> (the adult phase 3 trial) and the natural history from human challenge studies. The remaining
> literature was used to justify structure or to compare results against, not to fit parameters.


## 1. Within-host viral dynamics and modelling

The skeleton of this model — the target-cell-limited replication loop, the eclipse phase, the interferon-driven transition to refractory — follows the structure of the papers below as written. That β and p cannot be identified separately from titre data alone, and that only the product β·T₀·p is fixed by the growth rate, is also a point these papers make in common.

- Baccam P et al. Kinetics of influenza A virus infection in humans. *J Virol* 2006;80:7590-9. [PMID 16840338](https://pubmed.ncbi.nlm.nih.gov/16840338/)
- Smith AM et al. Influenza A virus infection kinetics: quantitative data and models. *Wiley Interdiscip Rev Syst Biol Med* 2011;3:429-45. [PMID 21197654](https://pubmed.ncbi.nlm.nih.gov/21197654/)
- Beauchemin CA et al. A review of mathematical models of influenza A infections within a host or cell culture: lessons learned and challenges ahead. *BMC Public Health* 2011;11 Suppl 1:S7. [PMID 21356136](https://pubmed.ncbi.nlm.nih.gov/21356136/)
- Canini L et al. Viral kinetic modeling: state of the art. *J Pharmacokinet Pharmacodyn* 2014;41:431-43. [PMID 24961742](https://pubmed.ncbi.nlm.nih.gov/24961742/)
- Canini L et al. Impact of different oseltamivir regimens on treating influenza A virus infection and resistance emergence: insights from a modelling study. *PLoS Comput Biol* 2014;10:e1003568. [PMID 24743564](https://pubmed.ncbi.nlm.nih.gov/24743564/)
- Sachak-Patwa R et al. A target-cell limited model can reproduce influenza infection dynamics in hosts with differing immune responses. *J Theor Biol* 2023;567:111491. [PMID 37044357](https://pubmed.ncbi.nlm.nih.gov/37044357/)
- Li K et al. Modelling the Effect of MUC1 on Influenza Virus Infection Kinetics and Macrophage Dynamics. *Viruses* 2021;13. [PMID 34066999](https://pubmed.ncbi.nlm.nih.gov/34066999/)
- Handel A et al. Neuraminidase inhibitor resistance in influenza: assessing the danger of its generation and spread. *PLoS Comput Biol* 2007;3:e240. [PMID 18069885](https://pubmed.ncbi.nlm.nih.gov/18069885/)

## 2. Human challenge studies and natural history

The calibration targets of A0 (peak titre · time of peak · duration of shedding · time of symptom onset · peak IL-6) all came from the papers in this section. Kaiser 2001 and Gentile 2001 are the basis of this model's SYM equation, in which symptoms follow IL-6 rather than virus.

- Carrat F et al. Time lines of infection and disease in human influenza: a review of volunteer challenge studies. *Am J Epidemiol* 2008;167:775-85. [PMID 18230677](https://pubmed.ncbi.nlm.nih.gov/18230677/)
- Hayden FG et al. Local and systemic cytokine responses during experimental human influenza A virus infection. Relation to symptom formation and host defense. *J Clin Invest* 1998;101:643-9. [PMID 9449698](https://pubmed.ncbi.nlm.nih.gov/9449698/)
- Fritz RS et al. Nasal cytokine and chemokine responses in experimental influenza A virus infection: results of a placebo-controlled trial of intravenous zanamivir treatment. *J Infect Dis* 1999;180:586-93. [PMID 10438343](https://pubmed.ncbi.nlm.nih.gov/10438343/)
- Hayden FG et al. Use of the oral neuraminidase inhibitor oseltamivir in experimental human influenza: randomized controlled trials for prevention and treatment. *JAMA* 1999;282:1240-6. [PMID 10517426](https://pubmed.ncbi.nlm.nih.gov/10517426/)
- Calfee DP et al. Safety and efficacy of once daily intranasal zanamivir in preventing experimental human influenza A infection. *Antivir Ther* 1999;4:143-9. [PMID 12731753](https://pubmed.ncbi.nlm.nih.gov/12731753/)
- Kaiser L et al. Symptom pathogenesis during acute influenza: interleukin-6 and other cytokine responses. *J Med Virol* 2001;64:262-8. [PMID 11424113](https://pubmed.ncbi.nlm.nih.gov/11424113/)
- Gentile DA et al. Effect of intranasal challenge with interleukin-6 on upper airway symptomatology and physiology in allergic and nonallergic patients. *Ann Allergy Asthma Immunol* 2001;86:531-6. [PMID 11379804](https://pubmed.ncbi.nlm.nih.gov/11379804/)

## 3. Baloxavir trials and pharmacology

The comparison figures in the A3 ledger (TTAS 80.2 / 53.8 / 53.7 h, cessation of shedding 96 / 72 / 24 h, day-2 titre −1.3 / −2.8 / −4.8 log₁₀) came from CAPSTONE-1. The argument for prophylactic dosing in A10 is BLOCKSTONE.

- Hayden FG et al. Baloxavir Marboxil for Uncomplicated Influenza in Adults and Adolescents. *N Engl J Med* 2018;379:913-923. [PMID 30184455](https://pubmed.ncbi.nlm.nih.gov/30184455/)
- Ison MG et al. Early treatment with baloxavir marboxil in high-risk adolescent and adult outpatients with uncomplicated influenza (CAPSTONE-2): a randomised, placebo-controlled, phase 3 trial. *Lancet Infect Dis* 2020;20:1204-1214. [PMID 32526195](https://pubmed.ncbi.nlm.nih.gov/32526195/)
- Ikematsu H et al. Baloxavir Marboxil for Prophylaxis against Influenza in Household Contacts. *N Engl J Med* 2020;383:309-320. [PMID 32640124](https://pubmed.ncbi.nlm.nih.gov/32640124/)
- Baker J et al. Baloxavir Marboxil Single-dose Treatment in Influenza-infected Children: A Randomized, Double-blind, Active Controlled Phase 3 Safety and Efficacy Trial (miniSTONE-2). *Pediatr Infect Dis J* 2020;39:700-705. [PMID 32516282](https://pubmed.ncbi.nlm.nih.gov/32516282/)
- Baker JB et al. Safety and Efficacy of Baloxavir Marboxil in Influenza-infected Children 5-11 Years of Age: A Post Hoc Analysis of a Phase 3 Study. *Pediatr Infect Dis J* 2023;42:983-989. [PMID 37595103](https://pubmed.ncbi.nlm.nih.gov/37595103/)
- Palmu S et al. A Phase 3 Safety and Efficacy Study of Baloxavir Marboxil in Children Less Than 1 Year Old With Suspected or Confirmed Influenza. *Pediatr Infect Dis J* 2025;44:645-649. [PMID 40279637](https://pubmed.ncbi.nlm.nih.gov/40279637/)
- Noshi T et al. In vitro characterization of baloxavir acid, a first-in-class cap-dependent endonuclease inhibitor of the influenza virus polymerase PA subunit. *Antiviral Res* 2018;160:109-117. [PMID 30316915](https://pubmed.ncbi.nlm.nih.gov/30316915/)
- Liu H et al. Simultaneous quantification of baloxavir marboxil and its active metabolite in human plasma using UHPLC-MS/MS: Application to a human pharmacokinetic study with different anticoagulants. *J Pharm Biomed Anal* 2024;249:116387. [PMID 39083919](https://pubmed.ncbi.nlm.nih.gov/39083919/)
- Ikematsu H et al. Comparative Effectiveness of Baloxavir Marboxil and Oseltamivir Treatment in Reducing Household Transmission of Influenza: A Post Hoc Analysis of the BLOCKSTONE Trial. *Influenza Other Respir Viruses* 2024;18:e13302. [PMID 38706384](https://pubmed.ncbi.nlm.nih.gov/38706384/)
- Ishiguro N et al. Clinical and virological outcomes with baloxavir compared with oseltamivir in pediatric patients aged 6 to < 12 years with influenza: an open-label, randomized, active-controlled trial protocol. *BMC Infect Dis* 2021;21:777. [PMID 34372769](https://pubmed.ncbi.nlm.nih.gov/34372769/)
- Ishiguro N et al. Clinical and Virologic Outcomes of Baloxavir Compared with Oseltamivir in Pediatric Patients with Influenza in Japan. *Infect Dis Ther* 2025;14:833-846. [PMID 40155497](https://pubmed.ncbi.nlm.nih.gov/40155497/)

## 4. Baloxavir resistance (PA/I38T)

The 50-fold shift in EC50 and the 18% replicative fitness cost used in the competitive release analysis of A6 came from this section. CAPSTONE-1 reported the emergence of PA/I38T · M · F in 9.7% of adults given baloxavir.

- Omoto S et al. Characterization of influenza virus variants induced by treatment with the endonuclease inhibitor baloxavir marboxil. *Sci Rep* 2018;8:9633. [PMID 29941893](https://pubmed.ncbi.nlm.nih.gov/29941893/)
- Uehara T et al. Treatment-Emergent Influenza Variant Viruses With Reduced Baloxavir Susceptibility: Impact on Clinical and Virologic Outcomes in Uncomplicated Influenza. *J Infect Dis* 2020;221:346-355. [PMID 31309975](https://pubmed.ncbi.nlm.nih.gov/31309975/)
- Lee LY et al. Evaluating the fitness of PA/I38T-substituted influenza A viruses with reduced baloxavir susceptibility in a competitive mixtures ferret model. *PLoS Pathog* 2021;17:e1009527. [PMID 33956888](https://pubmed.ncbi.nlm.nih.gov/33956888/)
- Kuroda T et al. In Vivo Antiviral Activity of Baloxavir against PA/I38T-Substituted Influenza A Viruses at Clinically Relevant Doses. *Viruses* 2023;15. [PMID 37243240](https://pubmed.ncbi.nlm.nih.gov/37243240/)
- Case JR et al. Impaired host shutoff is a fitness cost associated with baloxavir marboxil resistance mutations in influenza A virus PA/PA-X nuclease domain. *PLoS Pathog* 2026;22:e1013550. [PMID 41662358](https://pubmed.ncbi.nlm.nih.gov/41662358/)
- Sun B et al. Emergence of the novel PA-D27G mutation conferring reduced baloxavir susceptibility in influenza A viruses circulating in China, 2018-2025. *Emerg Microbes Infect* 2026;15:2620222. [PMID 41555529](https://pubmed.ncbi.nlm.nih.gov/41555529/)

## 5. Neuraminidase inhibitors

Both the PK of oseltamivir carboxylate (CL/F 18.8 L/h, ELF:plasma ≈ 1) and the dispute over the size of the clinical effect (Cochrane versus the Dobson IPD meta-analysis) are in this section. This model takes neither side and writes both values side by side in the ledger.

- Treanor JJ et al. Efficacy and safety of the oral neuraminidase inhibitor oseltamivir in treating acute influenza: a randomized controlled trial. US Oral Neuraminidase Study Group. *JAMA* 2000;283:1016-24. [PMID 10697061](https://pubmed.ncbi.nlm.nih.gov/10697061/)
- Nicholson KG et al. Efficacy and safety of oseltamivir in treatment of acute influenza: a randomised controlled trial. Neuraminidase Inhibitor Flu Treatment Investigator Group. *Lancet* 2000;355:1845-50. [PMID 10866439](https://pubmed.ncbi.nlm.nih.gov/10866439/)
- Dobson J et al. Oseltamivir treatment for influenza in adults: a meta-analysis of randomised controlled trials. *Lancet* 2015;385:1729-1737. [PMID 25640810](https://pubmed.ncbi.nlm.nih.gov/25640810/)
- Jefferson T et al. Oseltamivir for influenza in adults and children: systematic review of clinical study reports and summary of regulatory comments. *BMJ* 2014;348:g2545. [PMID 24811411](https://pubmed.ncbi.nlm.nih.gov/24811411/)
- Kaiser L et al. Impact of oseltamivir treatment on influenza-related lower respiratory tract complications and hospitalizations. *Arch Intern Med* 2003;163:1667-72. [PMID 12885681](https://pubmed.ncbi.nlm.nih.gov/12885681/)
- Kamal MA et al. Identification of new oral dosing regimens for the neuraminidase inhibitor oseltamivir in patients with moderate and severe renal impairment. *Clin Pharmacol Drug Dev* 2015;4:326-36. [PMID 27137141](https://pubmed.ncbi.nlm.nih.gov/27137141/)
- Su CP et al. Inhaled Zanamivir vs Oral Oseltamivir to Prevent Influenza-related Hospitalization or Death: A Nationwide Population-based Quasi-experimental Study. *Clin Infect Dis* 2022;75:1273-1279. [PMID 35299245](https://pubmed.ncbi.nlm.nih.gov/35299245/)
- Saisho Y et al. Pharmacokinetics and safety of intravenous peramivir, neuraminidase inhibitor of influenza virus, in healthy Japanese subjects. *Antivir Ther* 2017;22:313-323. [PMID 27805571](https://pubmed.ncbi.nlm.nih.gov/27805571/)
- Dong L et al. In Vitro and In Vivo Assessment of Pharmacokinetic Profile of Peramivir in the Context of Inhalation Therapy. *Pharmaceuticals (Basel)* 2025;18. [PMID 40005995](https://pubmed.ncbi.nlm.nih.gov/40005995/)
- Li L et al. Comparison of double-dose vs standard-dose oseltamivir in the treatment of influenza: A systematic review and meta-analysis. *J Clin Pharm Ther* 2020;45:918-926. [PMID 32497319](https://pubmed.ncbi.nlm.nih.gov/32497319/)
- Peng Y et al. Bioequivalence and Safety of Two Oseltamivir Phosphate for Oral Suspension in Healthy Chinese Subjects Under Fasting and Fed Conditions: A Randomized, Open‑Label, Single‑Dose, Crossover Study. *Clin Pharmacol Drug Dev* 2026;15:e1612. [PMID 41045035](https://pubmed.ncbi.nlm.nih.gov/41045035/)

## 6. NAI resistance and surveillance

The basis of the H275Y profile (EC50 300-fold, fitness cost 28%); Baz 2009, a case of resistance emerging during prophylaxis, bears directly on the A6 prediction that 'early dosing buys clinical benefit and selection for resistance at the same time'.

- Baz M et al. Emergence of oseltamivir-resistant pandemic H1N1 virus during prophylaxis. *N Engl J Med* 2009;361:2296-7. [PMID 19907034](https://pubmed.ncbi.nlm.nih.gov/19907034/)
- Takashita E et al. Global update on the susceptibilities of human influenza viruses to neuraminidase inhibitors and the cap-dependent endonuclease inhibitor baloxavir, 2017-2018. *Antiviral Res* 2020;175:104718. [PMID 32004620](https://pubmed.ncbi.nlm.nih.gov/32004620/)
- Farrukee R et al. Characterization of Influenza B Virus Variants with Reduced Neuraminidase Inhibitor Susceptibility. *Antimicrob Agents Chemother* 2018;62. [PMID 30201817](https://pubmed.ncbi.nlm.nih.gov/30201817/)
- Treurnicht FK et al. Replacement of neuraminidase inhibitor-susceptible influenza A(H1N1) with resistant phenotype in 2008 and circulation of susceptible influenza A and B viruses during 2009-2013, South Africa. *Influenza Other Respir Viruses* 2019;13:54-63. [PMID 30218485](https://pubmed.ncbi.nlm.nih.gov/30218485/)
- Duwe SC et al. Increase of Synergistic Secondary Antiviral Mutations in the Evolution of A(H1N1)pdm09 Influenza Virus Neuraminidases. *Viruses* 2024;16. [PMID 39066271](https://pubmed.ncbi.nlm.nih.gov/39066271/)
- Oh DY et al. Preparing for the Next Influenza Season: Monitoring the Emergence and Spread of Antiviral Resistance. *Infect Drug Resist* 2023;16:949-959. [PMID 36814825](https://pubmed.ncbi.nlm.nih.gov/36814825/)
- Ait-Aissa A et al. Surveillance for antiviral resistance among influenza viruses circulating in Algeria during five consecutive influenza seasons (2009-2014). *J Med Virol* 2018;90:844-853. [PMID 29315673](https://pubmed.ncbi.nlm.nih.gov/29315673/)
- Ngiam JN et al. Antiviral Resistance in Influenza: Clinical and Public Health Implications. *Influenza Other Respir Viruses* 2026;20:e70269. [PMID 42306964](https://pubmed.ncbi.nlm.nih.gov/42306964/)
- Cochin M et al. Characterization of oseltamivir-resistant A(H5N1) clade 2.3.4.4b, genotype D1.1 variants identified in poultry farms of British Columbia, Canada. *Emerg Microbes Infect* 2026;15:2686474. [PMID 42419257](https://pubmed.ncbi.nlm.nih.gov/42419257/)

## 7. Innate immunity and interferon

The basis of the T → R (refractory) transition, the IFN antagonism of NS1, and of 'suppression of antibacterial defence by interferon', one of the two conditions for secondary bacterial infection in A11.

- Shahangian A et al. Type I IFNs mediate development of postinfluenza bacterial pneumonia in mice. *J Clin Invest* 2009;119:1910-20. [PMID 19487810](https://pubmed.ncbi.nlm.nih.gov/19487810/)
- Cheung PH et al. Virus subtype-specific suppression of MAVS aggregation and activation by PB1-F2 protein of influenza A (H7N9) virus. *PLoS Pathog* 2020;16:e1008611. [PMID 32511263](https://pubmed.ncbi.nlm.nih.gov/32511263/)
- Ihazmade H et al. Impact of IFITM3 rs12252 and IFNAR2 rs2236757 SNP polymorphisms on influenza virus infection outcomes in Moroccan patients. *Mol Biol Rep* 2025;53:2. [PMID 41128776](https://pubmed.ncbi.nlm.nih.gov/41128776/)
- An S et al. Initial Influenza Virus Replication Can Be Limited in Allergic Asthma Through Rapid Induction of Type III Interferons in Respiratory Epithelium. *Front Immunol* 2018;9:986. [PMID 29867963](https://pubmed.ncbi.nlm.nih.gov/29867963/)
- Blümke J et al. Interferon lambda in anti-viral defense and cancer: dual roles, mechanism and therapeutic potential. *J Transl Med* 2026;24. [PMID 42092929](https://pubmed.ncbi.nlm.nih.gov/42092929/)

## 8. Adaptive immunity

The basis of the timing and magnitude of CD8 expansion (peak d6–8, 20–100-fold over naive), blockade of entry by mucosal IgA, and antibody-enhanced clearance of virions.

- Saidu A et al. Highly focused human CD8+ T-cell response in the lower airways during acute influenza infection. *J Immunol* 2026;215. [PMID 42152613](https://pubmed.ncbi.nlm.nih.gov/42152613/)
- Hertoghs N et al. A group 1 hemagglutinin stem vaccine elicits broad humoral responses against influenza in phase 1/2a study. *Nat Commun* 2026;17. [PMID 41974675](https://pubmed.ncbi.nlm.nih.gov/41974675/)
- Vanderven HA et al. Fc-mediated functions and the treatment of severe respiratory viral infections with passive immunotherapy - a balancing act. *Front Immunol* 2023;14:1307398. [PMID 38077353](https://pubmed.ncbi.nlm.nih.gov/38077353/)
- Papargyris L et al. Innate immune responsiveness predicts enhanced cellular immunity and symptomatic disease after controlled human influenza infection. *Nat Med* 2026;32:2556-2569. [PMID 42387215](https://pubmed.ncbi.nlm.nih.gov/42387215/)

## 9. Secondary bacterial infection

The BAC compartment of A11 is built as the product of two conditions — adhesion sites exposed by sloughed epithelium, and suppression of antibacterial defence by type I interferon. The threshold itself is an uncalibrated assumption, and the README says so.

- McCullers JA Insights into the interaction between influenza virus and pneumococcus. *Clin Microbiol Rev* 2006;19:571-82. [PMID 16847087](https://pubmed.ncbi.nlm.nih.gov/16847087/)
- Shahangian A et al. Type I IFNs mediate development of postinfluenza bacterial pneumonia in mice. *J Clin Invest* 2009;119:1910-20. [PMID 19487810](https://pubmed.ncbi.nlm.nih.gov/19487810/)
- Kato K et al. Influenza A virus infection induces initial proliferation of commensal Streptococcus pneumoniae in the larynx leading to dissemination into the lower respiratory tract. *J Virol* 2026;100:e0055526. [PMID 42370676](https://pubmed.ncbi.nlm.nih.gov/42370676/)
- Seo JH et al. Current insights into bacterial secondary infection following influenza A virus infection. *Front Microbiol* 2026;17:1851115. [PMID 42267108](https://pubmed.ncbi.nlm.nih.gov/42267108/)
- Malainou C et al. TNF superfamily member 14 drives post-influenza depletion of alveolar macrophages, enabling secondary pneumococcal pneumonia. *J Clin Invest* 2026;136. [PMID 41252214](https://pubmed.ncbi.nlm.nih.gov/41252214/)
- Alshammari AK et al. Understanding the Molecular Interactions Between Influenza A Virus and Streptococcus Proteins in Co-Infection: A Scoping Review. *Pathogens* 2025;14. [PMID 40005491](https://pubmed.ncbi.nlm.nih.gov/40005491/)
- Zhang H Concerns of using sialidase fusion protein as an experimental drug to combat seasonal and pandemic influenza. *J Antimicrob Chemother* 2008;62:219-23. [PMID 18238888](https://pubmed.ncbi.nlm.nih.gov/18238888/)

## 10. Severe influenza · the immunocompromised host · adjunctive therapy

The basis of A9 (corticosteroids) and A12 (the immunocompromised host). The observational mortality signal for steroids is reproduced in direction by this model from mechanism alone, with no mortality term.

- Lansbury L et al. Corticosteroids as adjunctive therapy in the treatment of influenza. *Cochrane Database Syst Rev* 2019;2:CD010406. [PMID 30798570](https://pubmed.ncbi.nlm.nih.gov/30798570/)
- Lansbury LE et al. Corticosteroids as Adjunctive Therapy in the Treatment of Influenza: An Updated Cochrane Systematic Review and Meta-analysis. *Crit Care Med* 2020;48:e98-e106. [PMID 31939808](https://pubmed.ncbi.nlm.nih.gov/31939808/)
- Muthuri SG et al. Effectiveness of neuraminidase inhibitors in reducing mortality in patients admitted to hospital with influenza A H1N1pdm09 virus infection: a meta-analysis of individual participant data. *Lancet Respir Med* 2014;2:395-404. [PMID 24815805](https://pubmed.ncbi.nlm.nih.gov/24815805/)
- Heldman MR et al. Influenza Antivirals for Prevention and Treatment in Immunocompromised People. *J Infect Dis* 2025;232:S243-S253. [PMID 41102612](https://pubmed.ncbi.nlm.nih.gov/41102612/)
- Hui DSC et al. Host Immunomodulatory Interventions in Severe Influenza. *J Infect Dis* 2025;232:S262-S272. [PMID 41102617](https://pubmed.ncbi.nlm.nih.gov/41102617/)
- Salvatore M et al. Baloxavir for the treatment of Influenza in allogeneic hematopoietic stem cell transplant recipients previously treated with oseltamivir. *Transpl Infect Dis* 2020;22:e13336. [PMID 32449254](https://pubmed.ncbi.nlm.nih.gov/32449254/)
- Euzen V et al. Zanamivir and baloxavir combination to cure persistent influenza and coronavirus infections after hematopoietic stem cell transplant. *Int J Antimicrob Agents* 2024;64:107281. [PMID 39047913](https://pubmed.ncbi.nlm.nih.gov/39047913/)
- Vetter P et al. Use of baloxavir as adjunctive antiviral therapy to neuraminidase inhibitors in severely immunocompromised individuals infected with influenza. *Antimicrob Agents Chemother* 2026;70:e0165925. [PMID 41891860](https://pubmed.ncbi.nlm.nih.gov/41891860/)
- Margaroli C et al. Spatial mapping of SARS-CoV-2 and H1N1 lung injury identifies differential transcriptional signatures. *Cell Rep Med* 2021;2:100242. [PMID 33778787](https://pubmed.ncbi.nlm.nih.gov/33778787/)

## 11. Other antivirals and biologics

The basis of favipiravir (lethal mutagenesis), peramivir, the M2 inhibitors (effectively abandoned because of S31N) and the anti-HA monoclonals. Of these, only oseltamivir and baloxavir are calibrated to clinical data; the rest are assumed values.

- Stevaert A et al. Nucleoside analogs for management of respiratory virus infections: mechanism of action and clinical efficacy. *Curr Opin Virol* 2022;57:101279. [PMID 36403338](https://pubmed.ncbi.nlm.nih.gov/36403338/)
- Geraghty RJ et al. Broad-Spectrum Antiviral Strategies and Nucleoside Analogues. *Viruses* 2021;13. [PMID 33924302](https://pubmed.ncbi.nlm.nih.gov/33924302/)
- Gregor J et al. Structural and Thermodynamic Analysis of the Resistance Development to Pimodivir (VX-787), the Clinical Inhibitor of Cap Binding to PB2 Subunit of Influenza A Polymerase. *Molecules* 2021;26. [PMID 33673017](https://pubmed.ncbi.nlm.nih.gov/33673017/)
- Stampolaki M et al. Adamantane-based inhibitors of the influenza A M2 proton channel: structure-based design, biological evaluation, and synthetic approaches. *RSC Med Chem* 2026;17:1811-1846. [PMID 41908668](https://pubmed.ncbi.nlm.nih.gov/41908668/)
- Dong J et al. Discovery of adamantane-based α-hydroxycarboxylic acid derivatives as potent M2-S31N blockers of influenza A virus. *Bioorg Chem* 2025;166:109124. [PMID 41151327](https://pubmed.ncbi.nlm.nih.gov/41151327/)
- Lim JJ et al. A Phase 2 Randomized, Double-Blind, Placebo-Controlled Trial of the Monoclonal Antibody MHAA4549A in Patients With Acute Uncomplicated Influenza A Infection. *Open Forum Infect Dis* 2022;9:ofab630. [PMID 35106315](https://pubmed.ncbi.nlm.nih.gov/35106315/)
- Mota KG et al. Monoclonal antibodies against influenza viruses: a clinical trials review. *Front Immunol* 2025;16:1669073. [PMID 41142796](https://pubmed.ncbi.nlm.nih.gov/41142796/)
- Wei J et al. Inhibition of cap-dependent endonuclease in influenza virus with ADC189: a pre-clinical analysis and phase I trial. *Front Med* 2025;19:347-358. [PMID 39832023](https://pubmed.ncbi.nlm.nih.gov/39832023/)

## 12. Endpoints · epidemiology · epithelial regeneration

The basis of the definition of TTAS (all seven symptoms mild or less for 21.5 hours), of the instruments that measure it, and of epithelial regeneration (LREG).

- Keeley TJH et al. Content validity and psychometric properties of the inFLUenza Patient-Reported Outcome Plus (FLU-PRO Plus(©)) instrument in patients with COVID-19. *Qual Life Res* 2023;32:1645-1657. [PMID 36703019](https://pubmed.ncbi.nlm.nih.gov/36703019/)
- Saito R et al. Duration of fever and symptoms in children after treatment with baloxavir marboxil and oseltamivir during the 2018-2019 season and detection of variant influenza a viruses with polymerase acidic subunit substitutions. *Antiviral Res* 2020;183:104951. [PMID 32987032](https://pubmed.ncbi.nlm.nih.gov/32987032/)
- Chung JR et al. Influenza vaccine effectiveness against outpatient acute respiratory illness with laboratory-confirmed influenza, United States, 2024-25 season. *Clin Infect Dis* 2026. [PMID 42442752](https://pubmed.ncbi.nlm.nih.gov/42442752/)
- Okoli GN et al. A systematic meta-analytic comparative evaluation of seasonal influenza vaccine effectiveness from test-negative design studies in the Northern Hemisphere pre/post COVID-19 pandemic. *Vaccine* 2026;88:128956. [PMID 42486051](https://pubmed.ncbi.nlm.nih.gov/42486051/)
- Chan LYH et al. Understanding spatiotemporal clustering of seasonal influenza in the United States. *BMC Infect Dis* 2026;26. [PMID 41781891](https://pubmed.ncbi.nlm.nih.gov/41781891/)
- Brcko IC et al. Phylodynamic reconstruction of H1N1pdm09 influenza virus transmission in Brazil: a decade of evolutionary dynamics. *Emerg Microbes Infect* 2026;15:2620237. [PMID 41555524](https://pubmed.ncbi.nlm.nih.gov/41555524/)
- Grieco L et al. Exploring the role of mass immunisation in influenza pandemic preparedness: A modelling study for the UK context. *Vaccine* 2020;38:5163-5170. [PMID 32576461](https://pubmed.ncbi.nlm.nih.gov/32576461/)
- Gagneux P et al. Human-specific regulation of alpha 2-6-linked sialic acids. *J Biol Chem* 2003;278:48245-50. [PMID 14500706](https://pubmed.ncbi.nlm.nih.gov/14500706/)
- Zeng H et al. Tropism and Infectivity of a Seasonal A(H1N1) and a Highly Pathogenic Avian A(H5N1) Influenza Virus in Primary Differentiated Ferret Nasal Epithelial Cell Cultures. *J Virol* 2019;93. [PMID 30814288](https://pubmed.ncbi.nlm.nih.gov/30814288/)
- Weiner AI et al. ΔNp63 drives dysplastic alveolar remodeling and restricts epithelial plasticity upon severe lung injury. *Cell Rep* 2022;41:111805. [PMID 36516758](https://pubmed.ncbi.nlm.nih.gov/36516758/)
- Lu T et al. Dysplastic epithelial repair promotes the tissue residence of lymphocytes to inhibit alveolar regeneration post viral infection. *Cell Stem Cell* 2026;33:108-124.e6. [PMID 41443194](https://pubmed.ncbi.nlm.nih.gov/41443194/)
- Yu Y et al. Timing-specific efficacy of antiviral baloxavir and anti-inflammatory oclacitinib monotherapies, and the benefits of their combination in treating influenza in mice. *Front Microbiol* 2026;17:1741128. [PMID 41918532](https://pubmed.ncbi.nlm.nih.gov/41918532/)


---

## What this model fails to reproduce (reported discrepancies)

Points of disagreement with the literature are not hidden but reported, numbered, inside the
analysis. All of them appear as they stand in the output of `flu_reference_check.py`.

| # | Where | Content |
|---|---|---|
| 1 | A2 | with a Hill slope of 1 (Michaelis) the 24-hour titre fall of CAPSTONE-1 (−4.8 log₁₀) cannot be produced **by any Emax**, because the lower bound on the residual production fraction is EC50/(C+EC50). The trial data are evidence about the **shape of the concentration-response curve**, not about potency. |
| 2 | A12 | the prediction that late antiviral treatment retains value in the immunocompromised host — consistent with clinical practice, but no randomised trial has tested it. It is merely unrefuted, not validated. |
| 3 | A8 | if symptoms track virus at all, baloxavir should give a greater symptomatic benefit than oseltamivir. CAPSTONE-1 could not separate them — 53.7 versus 53.8 hours — despite a 2 log difference in day-2 titre. |
| 4 | A6 | because the model is deterministic, a resistant lineage at the 10⁻⁴ level 'exists' in every virtual patient. Real I38T emergence is a stochastic establishment event occurring in about 10% of adults. |
| 5 | A7 | **the ordering comes out backwards.** In this model pre-existing immunity slows the wild type and so leaves more of the target-cell field, which increases selection for resistance instead. That is the opposite direction to the observation that I38T is commoner in children, and it is the clearest signal that the well-mixed epithelium assumption is wrong for this question. |

---

## Tools and methodology

- mrgsolve — <https://mrgsolve.org/>
- Graphviz — <https://graphviz.org/>
- PubMed E-utilities (used to verify this list) — <https://www.ncbi.nlm.nih.gov/books/NBK25501/>

**97 papers in total** (all confirmed by PubMed lookup).

> This is a model for education and research. It must not be used for clinical decision-making.
