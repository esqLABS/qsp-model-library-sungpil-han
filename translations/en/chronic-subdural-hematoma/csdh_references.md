# Chronic Subdural Haematoma (cSDH) QSP Model — References
# Chronic Subdural Haematoma — QSP Model References

A total of **130** references. Every PMID was checked programmatically against PubMed E-utilities to confirm it actually exists (`csdh_fetch_references.py`). To see which literature each of the model's structural claims relies on, see the section breakdown below together with the "What comes from the literature and what is the model's claim" section in `README.md`.

All PMIDs were retrieved programmatically via the PubMed E-utilities API and therefore correspond to real records. Numbers quoted in the model comments come from `csdh_verify.py`, not from these papers directly.

---

## Landmark trials the model is calibrated against

| Trial | Design | Value used in the model |
|------|------|--------------------|
| **Santarius 2009** (Lancet) | Randomised subdural drain after burr-hole drainage | recurrence 24% (no drain) vs 9.3% (drain) |
| **ATOCH 2018** (JAMA Neurol) | Atorvastatin 20 mg for 8 weeks, conservative treatment | conversion to surgery 11.2% vs 23.5%; neurological improvement 45.9% vs 28.6% |
| **Dex-CSDH 2020** (NEJM) | Dexamethasone 2-week tapering regimen, n=748 | favourable outcome 83.9% vs 90.3%; reoperation 1.7% vs 7.1%; SAE 16.0% vs 6.4% |
| **EMBOLISE 2024** (NEJM) | Surgery-adjunct MMA embolisation, 180 days | treatment failure 4.1% vs 11.3% |
| **MAGIC-MT 2024** | MMA embolisation, 90-day symptomatic recurrence | 6.7% vs 9.9% (p=0.10, negative); SAE 6.7% vs 11.6% |

> EMBOLISE (180 days, surgery-adjunct) is positive and MAGIC-MT (90 days) is negative. Rather than build that discrepancy in
> as an assumption, this model **derives it from structure**: because embolisation acts on the SOURCE rather than the STOCK,
> its effect is an integral over the neovascular decay time constant (~9–18 days), so at the 30–90-day time point it is
> still small and grows large by 180 days. See `CSDH_mmae_timecourse()` for the detailed calculation.

---

## Epidemiology & Natural History

- Tommiska P et al. (2022). *Incidence of surgery for chronic subdural hematoma in Finland during 1997-2014: a nationwide study*. J Neurosurg. [PMID 34507291](https://pubmed.ncbi.nlm.nih.gov/34507291/)
- Huang J et al. (2020). *Drug treatment of chronic subdural hematoma*. Expert Opin Pharmacother. [PMID 31957506](https://pubmed.ncbi.nlm.nih.gov/31957506/)
- Dziho A et al. (2025). *Global prevalence and incidence of chronic subdural hematoma: A systematic review*. Brain Spine. [PMID 41439175](https://pubmed.ncbi.nlm.nih.gov/41439175/)
- Stubbs DJ et al. (2021). *Incidence of chronic subdural haematoma: a single-centre exploration of the effects of an ageing population with a review of the literature*. Acta Neurochir (Wien). [PMID 34181085](https://pubmed.ncbi.nlm.nih.gov/34181085/)
- Soleman J et al. (2017). *The conservative and pharmacological management of chronic subdural haematoma*. Swiss Med Wkly. [PMID 28102879](https://pubmed.ncbi.nlm.nih.gov/28102879/)
- Liu LX et al. (2019). *Risk Factors for Recurrence of Chronic Subdural Hematoma: A Single Center Experience*. World Neurosurg. [PMID 31450003](https://pubmed.ncbi.nlm.nih.gov/31450003/)
- You W et al. (2018). *Prevalence of and risk factors for recurrence of chronic subdural hematoma*. Acta Neurochir (Wien). [PMID 29532258](https://pubmed.ncbi.nlm.nih.gov/29532258/)
- Hamou H et al. (2022). *Risk factors of recurrence in chronic subdural hematoma and a proposed extended classification of internal architecture as a predictor of recurrence*. Neurosurg Rev. [PMID 35461433](https://pubmed.ncbi.nlm.nih.gov/35461433/)
- Akyuz ME et al. (2023). *Identification of Radiological and Clinical Factors that Increase the Risk of Chronic Subdural Hematoma Recurrence*. Turk Neurosurg. [PMID 37885309](https://pubmed.ncbi.nlm.nih.gov/37885309/)

## Pathophysiology: Neomembrane & Neovascularisation

- Killeffer JA et al. (2000). *The outer neomembrane of chronic subdural hematoma*. Neurosurg Clin N Am. [PMID 10918009](https://pubmed.ncbi.nlm.nih.gov/10918009/)
- Nanko N et al. (2009). *Involvement of hypoxia-inducible factor-1alpha and vascular endothelial growth factor in the mechanism of development of chronic subdural hematoma*. Neurol Med Chir (Tokyo). [PMID 19779281](https://pubmed.ncbi.nlm.nih.gov/19779281/)
- Hohenstein A et al. (2005). *Increased mRNA expression of VEGF within the hematoma and imbalance of angiopoietin-1 and -2 mRNA within the neomembranes of chronic subdural hematoma*. J Neurotrauma. [PMID 15892598](https://pubmed.ncbi.nlm.nih.gov/15892598/)
- Lee KS et al. (1998). *Origin of chronic subdural haematoma and relation to traumatic subdural lesions*. Brain Inj. [PMID 9839025](https://pubmed.ncbi.nlm.nih.gov/9839025/)
- Takei J et al. (2021). *Significantly high concentrations of vascular endothelial growth factor in chronic subdural hematoma with trabecular formation*. Clin Neurol Neurosurg. [PMID 33545457](https://pubmed.ncbi.nlm.nih.gov/33545457/)
- Edlmann E et al. (2022). *Dexamethasone reduces vascular endothelial growth factor in comparison to placebo in post-operative chronic subdural hematoma samples: A target for future drug therapy?*. Front Neurol. [PMID 36158966](https://pubmed.ncbi.nlm.nih.gov/36158966/)
- Shono T et al. (2001). *Vascular endothelial growth factor in chronic subdural haematomas*. J Clin Neurosci. [PMID 11535006](https://pubmed.ncbi.nlm.nih.gov/11535006/)
- Li F et al. (2017). *Correlation of vascular endothelial growth factor with magnetic resonance imaging in chronic subdural hematomas*. J Neurol Sci. [PMID 28477686](https://pubmed.ncbi.nlm.nih.gov/28477686/)
- Weigel R et al. (2014). *Vascular endothelial growth factor concentration in chronic subdural hematoma fluid is related to computed tomography appearance and exudation rate*. J Neurotrauma. [PMID 24245657](https://pubmed.ncbi.nlm.nih.gov/24245657/)

## Pathophysiology: Inflammation & Haem

- Frati A et al. (2004). *Inflammation markers and risk factors for recurrence in 35 patients with a posttraumatic chronic subdural hematoma: a prospective study*. J Neurosurg. [PMID 14743908](https://pubmed.ncbi.nlm.nih.gov/14743908/)
- Yamashima T et al. (1984). *Clinicopathological study of acute subdural haematoma in the chronic healing stage. Clinical, histological and ultrastructural comparisons with chronic subdural haematoma*. Neurochirurgia (Stuttg). [PMID 6483072](https://pubmed.ncbi.nlm.nih.gov/6483072/)
- Hua C et al. (2016). *Role of Matrix Metalloproteinase-2, Matrix Metalloproteinase-9, and Vascular Endothelial Growth Factor in the Development of Chronic Subdural Hematoma*. J Neurotrauma. [PMID 25646653](https://pubmed.ncbi.nlm.nih.gov/25646653/)
- Trieu M et al. (2026). *Pathophysiology of chronic subdural hematoma-new insights*. Ther Adv Neurol Disord. [PMID 41953848](https://pubmed.ncbi.nlm.nih.gov/41953848/)
- Schack A et al. (2025). *Immunoprofile of Radiologic Chronic Subdural Hematoma Subtypes*. World Neurosurg. [PMID 39461416](https://pubmed.ncbi.nlm.nih.gov/39461416/)
- Wang ZF et al. (2010). *[Role of matrix metalloproteinase in transformation of subdural effusion into chronic subdural hematoma]*. Nan Fang Yi Ke Da Xue Xue Bao. [PMID 20501425](https://pubmed.ncbi.nlm.nih.gov/20501425/)

## Pathophysiology: Local Hyperfibrinolysis

- Osuka K et al. (2026). *Roles of tissue plasminogen activator in chronic subdural hematoma, independent of the fibrinolytic system*. Neurosurg Rev. [PMID 41634399](https://pubmed.ncbi.nlm.nih.gov/41634399/)
- Fujisawa H et al. (1991). *Immunohistochemical localization of tissue-type plasminogen activator in the lining wall of chronic subdural hematoma*. Surg Neurol. [PMID 1905064](https://pubmed.ncbi.nlm.nih.gov/1905064/)
- Ito H et al. (1976). *Role of local hyperfibrinolysis in the etiology of chronic subdural hematoma*. J Neurosurg. [PMID 132513](https://pubmed.ncbi.nlm.nih.gov/132513/)
- Ito H et al. (1978). *Fibrinolytic enzyme in the lining walls of chronic subdural hematoma*. J Neurosurg. [PMID 146730](https://pubmed.ncbi.nlm.nih.gov/146730/)
- Suzuki M et al. (1998). *Local hypercoagulative activity precedes hyperfibrinolytic activity in the subdural space during development of chronic subdural haematoma from subdural effusion*. Acta Neurochir (Wien). [PMID 9638263](https://pubmed.ncbi.nlm.nih.gov/9638263/)
- Cuny E (2001). *[Physiopathology of chronic subdural hematoma]*. Neurochirurgie. [PMID 11915758](https://pubmed.ncbi.nlm.nih.gov/11915758/)

## Fluid Dynamics & Clearance

- Labadie EL et al. (1974). *Chronic subdural hematoma: concepts of physiopathogenesis. A review*. Can J Neurol Sci. [PMID 4613449](https://pubmed.ncbi.nlm.nih.gov/4613449/)
- Jacob L et al. (2022). *Conserved meningeal lymphatic drainage circuits in mice and humans*. J Exp Med. [PMID 35776089](https://pubmed.ncbi.nlm.nih.gov/35776089/)
- Louveau A et al. (2017). *Understanding the functions and relationships of the glymphatic system and meningeal lymphatics*. J Clin Invest. [PMID 28862640](https://pubmed.ncbi.nlm.nih.gov/28862640/)
- McDonald DM et al. (2025). *Cerebrospinal fluid draining lymphatics in health and disease: advances and controversies*. Nat Cardiovasc Res. [PMID 40921861](https://pubmed.ncbi.nlm.nih.gov/40921861/)
- Lohat V et al. (2025). *Distinct Role of Dural and Leptomeningeal Macrophages in Maintaining Cerebrospinal Fluid Drainage to Meningeal Lymphatic Vessels*. Am J Pathol. [PMID 40541712](https://pubmed.ncbi.nlm.nih.gov/40541712/)
- Chen Y et al. (2024). *Vitamin D accelerates the subdural hematoma clearance through improving the meningeal lymphatic vessel function*. Mol Cell Biochem. [PMID 38294731](https://pubmed.ncbi.nlm.nih.gov/38294731/)
- Yuan J et al. (2024). *Inactivation of ERK1/2 signaling mediates dysfunction of basal meningeal lymphatic vessels in experimental subdural hematoma*. Theranostics. [PMID 38164141](https://pubmed.ncbi.nlm.nih.gov/38164141/)
- Liu X et al. (2020). *Subdural haematomas drain into the extracranial lymphatic system through the meningeal lymphatic vessels*. Acta Neuropathol Commun. [PMID 32059751](https://pubmed.ncbi.nlm.nih.gov/32059751/)

## Intracranial Mechanics & Brain Re-expansion

- Juan WS et al. (2012). *Multiple tenting techniques improve dead space obliteration in the surgical treatment for patients with giant calcified chronic subdural hematoma*. Acta Neurochir (Wien). [PMID 22109694](https://pubmed.ncbi.nlm.nih.gov/22109694/)
- Ishiwata Y et al. (1987). *[Subdural tension pneumocephalus following surgery of chronic subdural hematoma]*. No Shinkei Geka. [PMID 3614535](https://pubmed.ncbi.nlm.nih.gov/3614535/)
- Shen J et al. (2019). *Clinical and radiological factors predicting recurrence of chronic subdural hematoma: A retrospective cohort study*. Injury. [PMID 31445831](https://pubmed.ncbi.nlm.nih.gov/31445831/)
- Takei J et al. (2025). *Differences in neuroradiological impacts of hematoma volume and midline shift on clinical symptoms and recurrence rate in patients with unilateral chronic subdural hematoma*. J Clin Neurosci. [PMID 39986185](https://pubmed.ncbi.nlm.nih.gov/39986185/)
- Kujiraoka Y et al. (2004). *Shaken baby syndrome manifesting as chronic subdural hematoma: importance of single photon emission computed tomography for treatment indications--case report*. Neurol Med Chir (Tokyo). [PMID 15347212](https://pubmed.ncbi.nlm.nih.gov/15347212/)
- Khorasanizadeh M et al. (2024). *The effect of patient age on the degree of midline shift caused by chronic subdural hematomas: a volumetric analysis*. J Neurosurg. [PMID 37877977](https://pubmed.ncbi.nlm.nih.gov/37877977/)
- Brasil S et al. (2025). *Monro-Kellie 4.0: moving from intracranial pressure to intracranial dynamics*. Crit Care. [PMID 40474297](https://pubmed.ncbi.nlm.nih.gov/40474297/)
- Brasil S et al. (2026). *New concepts in intracranial compliance: pathophysiology, monitoring and clinical implications*. Curr Opin Crit Care. [PMID 41384504](https://pubmed.ncbi.nlm.nih.gov/41384504/)
- Kalisvaart ACJ et al. (2020). *An update to the Monro-Kellie doctrine to reflect tissue compliance after severe ischemic and hemorrhagic stroke*. Sci Rep. [PMID 33328490](https://pubmed.ncbi.nlm.nih.gov/33328490/)

## Surgical Management

- Wan Y et al. (2019). *Single Versus Double Burr Hole Craniostomy in Surgical Treatment of Chronic Subdural Hematoma: A Meta-Analysis*. World Neurosurg. [PMID 31323397](https://pubmed.ncbi.nlm.nih.gov/31323397/)
- Kutty SA et al. (2014). *Chronic subdural hematoma: a comparison of recurrence rates following burr-hole craniostomy with and without drains*. Turk Neurosurg. [PMID 25050672](https://pubmed.ncbi.nlm.nih.gov/25050672/)
- Wu L et al. (2023). *Exhaustive drainage versus fixed-time drainage for chronic subdural hematoma after one-burr hole craniostomy (ECHO): study protocol for a multicenter randomized controlled trial*. Trials. [PMID 36941714](https://pubmed.ncbi.nlm.nih.gov/36941714/)
- Mahmoodkhani M et al. (2022). *Half-Saline Versus Normal-Saline as Irrigation Solutions in Burr Hole Craniostomy to Treat Chronic Subdural Hematomata: A Randomized Clinical Trial*. Korean J Neurotrauma. [PMID 36381457](https://pubmed.ncbi.nlm.nih.gov/36381457/)
- Hutchinson PJ et al. (2024). *A randomised, double blind, placebo-controlled trial of a two-week course of dexamethasone for adult patients with a symptomatic Chronic Subdural Haematoma (Dex-CSDH trial)*. Health Technol Assess. [PMID 38512045](https://pubmed.ncbi.nlm.nih.gov/38512045/)
- Santarius T et al. (2009). *Use of drains versus no drains after burr-hole evacuation of chronic subdural haematoma: a randomised controlled trial*. Lancet. [PMID 19782872](https://pubmed.ncbi.nlm.nih.gov/19782872/)
- Guilfoyle MR et al. (2017). *Improved long-term survival with subdural drains following evacuation of chronic subdural haematoma*. Acta Neurochir (Wien). [PMID 28349381](https://pubmed.ncbi.nlm.nih.gov/28349381/)
- Aljabali A et al. (2024). *Irrigation versus no irrigation in the treatment of chronic subdural hematoma: An updated systematic review and meta-analysis of 1581 patients*. Neurosurg Rev. [PMID 38538863](https://pubmed.ncbi.nlm.nih.gov/38538863/)
- Raj R et al. (2024). *Burr-hole drainage with or without irrigation for chronic subdural haematoma (FINISH): a Finnish, nationwide, parallel-group, multicentre, randomised, controlled, non-inferiority trial*. Lancet. [PMID 38852600](https://pubmed.ncbi.nlm.nih.gov/38852600/)
- Schack A et al. (2025). *Irrigation practices in surgical evacuation of chronic subdural hematoma: systematic review and meta-analysis of technique, fluid type, and temperature*. Neurosurg Rev. [PMID 41217544](https://pubmed.ncbi.nlm.nih.gov/41217544/)
- Hellwig D et al. (1996). *Endoscopic treatment of septated chronic subdural hematoma*. Surg Neurol. [PMID 8638225](https://pubmed.ncbi.nlm.nih.gov/8638225/)
- Shahid AH et al. (2025). *Endoscopic subdural membranectomy for multi-septated chronic subdural hematoma: Finding a safe solution when middle meningeal artery embolization is not feasible*. Surg Neurol Int. [PMID 40469338](https://pubmed.ncbi.nlm.nih.gov/40469338/)
- Khaleghi M et al. (2026). *The role of endoscope-assisted septectomy and membranectomy for complex chronic subdural hematomas: safety, efficacy, technical feasibility. Patient series*. J Neurosurg Case Lessons. [PMID 41871398](https://pubmed.ncbi.nlm.nih.gov/41871398/)

## Middle Meningeal Artery Embolisation

- Davies JM et al. (2024). *Adjunctive Middle Meningeal Artery Embolization for Subdural Hematoma*. N Engl J Med. [PMID 39565988](https://pubmed.ncbi.nlm.nih.gov/39565988/)
- Fiorella D et al. (2025). *Embolization of the Middle Meningeal Artery for Chronic Subdural Hematoma*. N Engl J Med. [PMID 39565980](https://pubmed.ncbi.nlm.nih.gov/39565980/)
- Gillespie CS et al. (2025). *Middle meningeal artery embolization for chronic subdural hematoma: meta-analysis of three randomized controlled trials and review of ongoing trials*. Acta Neurochir (Wien). [PMID 40493076](https://pubmed.ncbi.nlm.nih.gov/40493076/)
- de Almeida Monteiro G et al. (2025). *Middle meningeal artery embolization for chronic subdural hematoma: a meta-analysis of randomized controlled trials with trial sequential analysis*. Neurosurg Rev. [PMID 40210750](https://pubmed.ncbi.nlm.nih.gov/40210750/)
- Liu J et al. (2024). *Middle Meningeal Artery Embolization for Nonacute Subdural Hematoma*. N Engl J Med. [PMID 39565989](https://pubmed.ncbi.nlm.nih.gov/39565989/)
- Shotar E et al. (2025). *Meningeal Embolization for Preventing Chronic Subdural Hematoma Recurrence After Surgery: The EMPROTECT Randomized Clinical Trial*. JAMA. [PMID 40471557](https://pubmed.ncbi.nlm.nih.gov/40471557/)
- Ironside N et al. (2021). *Middle meningeal artery embolization for chronic subdural hematoma: a systematic review and meta-analysis*. J Neurointerv Surg. [PMID 34193592](https://pubmed.ncbi.nlm.nih.gov/34193592/)
- Siddiq F et al. (2025). *Consensus Statement on Middle Meningeal Artery Embolization in Chronic Subdural Hematoma Treatment: A Guideline from the Society of Vascular and Interventional Neurology Guidelines and Practice Standards Committee*. Stroke Vasc Interv Neurol. [PMID 41608698](https://pubmed.ncbi.nlm.nih.gov/41608698/)
- Shankar J et al. (2025). *Embolization of middle meningeal artery for chronic subdural hematoma: Do we have sufficient evidence?*. Interv Neuroradiol. [PMID 38592031](https://pubmed.ncbi.nlm.nih.gov/38592031/)
- Chida K et al. (2026). *Clinical and Radiographic Course Following Middle Meningeal Artery Embolization for Chronic Subdural Hematoma*. J Neuroendovasc Ther. [PMID 41716851](https://pubmed.ncbi.nlm.nih.gov/41716851/)
- Jia Q et al. (2026). *Dangerous Ophthalmic Anastomoses During Middle Meningeal Artery Embolization: Anatomical Pitfalls, Hemodynamic Mechanisms, and an Illustrative Case of Delayed Visual Loss*. Clin Neuroradiol. [PMID 42243540](https://pubmed.ncbi.nlm.nih.gov/42243540/)
- Sattur MG et al. (2020). *Anomalous "Middle" Meningeal Artery from Basilar Artery and Implications for Neuroendovascular Surgery: Case Report and Review of Literature*. World Neurosurg. [PMID 31574331](https://pubmed.ncbi.nlm.nih.gov/31574331/)
- Manelfe C et al. (1986). *Preoperative embolization of intracranial meningiomas*. AJNR Am J Neuroradiol. [PMID 3096121](https://pubmed.ncbi.nlm.nih.gov/3096121/)

## Corticosteroids

- Hutchinson PJ et al. (2020). *Trial of Dexamethasone for Chronic Subdural Hematoma*. N Engl J Med. [PMID 33326713](https://pubmed.ncbi.nlm.nih.gov/33326713/)
- Miah IP et al. (2023). *Dexamethasone versus Surgery for Chronic Subdural Hematoma*. N Engl J Med. [PMID 37314705](https://pubmed.ncbi.nlm.nih.gov/37314705/)
- Edlmann E et al. (2019). *Dex-CSDH randomised, placebo-controlled trial of dexamethasone for chronic subdural haematoma: report of the internal pilot phase*. Sci Rep. [PMID 30971773](https://pubmed.ncbi.nlm.nih.gov/30971773/)
- Jiang RC et al. (2021). *Atorvastatin combined with dexamethasone in chronic subdural haematoma (ATOCH II): study protocol for a randomized controlled trial*. Trials. [PMID 34895306](https://pubmed.ncbi.nlm.nih.gov/34895306/)
- Delgado-López PD et al. (2009). *Dexamethasone treatment in chronic subdural haematoma*. Neurocirugia (Astur). [PMID 19688136](https://pubmed.ncbi.nlm.nih.gov/19688136/)
- Holl DC et al. (2019). *Corticosteroid treatment compared with surgery in chronic subdural hematoma: a systematic review and meta-analysis*. Acta Neurochir (Wien). [PMID 30972566](https://pubmed.ncbi.nlm.nih.gov/30972566/)
- Tang G et al. (2021). *The Efficacy of Adjuvant Corticosteroids in Surgical Management of Chronic Subdural Hematoma: A Systematic Review and Meta-Analysis*. Front Neurol. [PMID 35095713](https://pubmed.ncbi.nlm.nih.gov/35095713/)
- Almenawer SA et al. (2014). *Chronic subdural hematoma management: a systematic review and meta-analysis of 34,829 patients*. Ann Surg. [PMID 24096761](https://pubmed.ncbi.nlm.nih.gov/24096761/)
- Villar J et al. (2020). *Dexamethasone treatment for the acute respiratory distress syndrome: a multicentre, randomised controlled trial*. Lancet Respir Med. [PMID 32043986](https://pubmed.ncbi.nlm.nih.gov/32043986/)
- Bateman RM et al. (2016). *36th International Symposium on Intensive Care and Emergency Medicine : Brussels, Belgium. 15-18 March 2016*. Crit Care. [PMID 27885969](https://pubmed.ncbi.nlm.nih.gov/27885969/)
- Ferreri AJ et al. (2016). *Chemoimmunotherapy with methotrexate, cytarabine, thiotepa, and rituximab (MATRix regimen) in patients with primary CNS lymphoma: results of the first randomisation of the International Extranodal Lymphoma Study Group-32 (IELSG32) phase 2 trial*. Lancet Haematol. [PMID 27132696](https://pubmed.ncbi.nlm.nih.gov/27132696/)

## Atorvastatin

- Jiang R et al. (2018). *Safety and Efficacy of Atorvastatin for Chronic Subdural Hematoma in Chinese Patients: A Randomized ClinicalTrial*. JAMA Neurol. [PMID 30073290](https://pubmed.ncbi.nlm.nih.gov/30073290/)
- Monteiro GA et al. (2024). *Efficacy and Safety of Atorvastatin for Chronic Subdural Hematoma: An Updated Systematic Review and Meta-Analysis*. World Neurosurg. [PMID 38759787](https://pubmed.ncbi.nlm.nih.gov/38759787/)
- Jiang R et al. (2015). *Effect of ATorvastatin On Chronic subdural Hematoma (ATOCH): a study protocol for a randomized controlled trial*. Trials. [PMID 26581842](https://pubmed.ncbi.nlm.nih.gov/26581842/)
- Wang D et al. (2014). *Effects of atorvastatin on chronic subdural hematoma: a preliminary report from three medical centers*. J Neurol Sci. [PMID 24269089](https://pubmed.ncbi.nlm.nih.gov/24269089/)
- Tian Y et al. (2022). *Establishment and validation of a prediction model for self-absorption probability of chronic subdural hematoma*. Front Neurol. [PMID 35937067](https://pubmed.ncbi.nlm.nih.gov/35937067/)

## Tranexamic Acid & Antithrombotics

- Lodewijkx R et al. (2021). *Tranexamic acid for chronic subdural hematoma*. Br J Neurosurg. [PMID 34334070](https://pubmed.ncbi.nlm.nih.gov/34334070/)
- Messias BR et al. (2025). *Chronic Subdural Hematoma and Tranexamic Acid: A Systematic Review*. Turk Neurosurg. [PMID 40577507](https://pubmed.ncbi.nlm.nih.gov/40577507/)
- Pan W et al. (2024). *Effectiveness of tranexamic acid on chronic subdural hematoma recurrence: a meta-analysis and systematic review*. Front Neurol. [PMID 38711565](https://pubmed.ncbi.nlm.nih.gov/38711565/)
- Albalkhi I et al. (2024). *Adjuvant Tranexamic Acid for Reducing Postoperative Recurrence of Chronic Subdural Hematoma in the Elderly: A Systematic Review and Meta-Analysis*. World Neurosurg. [PMID 38101544](https://pubmed.ncbi.nlm.nih.gov/38101544/)
- Mishra R et al. (2025). *Role of Tranexamic Acid in the Management of Chronic Subdural Hematoma: A Systematic Review and Meta-analysis of Randomized Controlled Trials*. Neurol India. [PMID 40652463](https://pubmed.ncbi.nlm.nih.gov/40652463/)
- Ducruet AF et al. (2012). *The surgical management of chronic subdural hematoma*. Neurosurg Rev. [PMID 21909694](https://pubmed.ncbi.nlm.nih.gov/21909694/)
- Raj R et al. (2025). *Restarting anticoagulation early versus late in patients with chronic subdural hematoma and atrial fibrillation (RELACS): a phase III international multicenter, randomized controlled, two-arm, assessor-blinded trial*. Trials. [PMID 41250110](https://pubmed.ncbi.nlm.nih.gov/41250110/)
- Dakhel S et al. (2026). *Timing of antithrombotic therapy resumption after surgical evacuation of chronic subdural hematoma: a propensity score-matched analysis*. Neurosurg Rev. [PMID 41912884](https://pubmed.ncbi.nlm.nih.gov/41912884/)
- Harland TA et al. (2026). *Treatment Durability After Subdural Evacuating Port System, Middle Meningeal Artery Embolization, and Combined Therapy for Chronic Subdural Hematoma*. Neurosurgery. [PMID 42240333](https://pubmed.ncbi.nlm.nih.gov/42240333/)
- Yadav YR et al. (2016). *Chronic subdural hematoma*. Asian J Neurosurg. [PMID 27695533](https://pubmed.ncbi.nlm.nih.gov/27695533/)
- Lee SW et al. (2024). *Risk Factors for the Recurrence of Chronic Subdural Hematoma*. Korean J Neurotrauma. [PMID 39021754](https://pubmed.ncbi.nlm.nih.gov/39021754/)
- SeyedAlinaghi S et al. (2025). *Risk Factors of Chronic Subdural Hematoma Recurrence: A Systematic Review*. Med J Islam Repub Iran. [PMID 41194951](https://pubmed.ncbi.nlm.nih.gov/41194951/)

## Pharmacokinetics

- Hoffmann M et al. (2025). *Elranatamab for Relapsed/Refractory Multiple Myeloma With Severe Renal Impairment Requiring Hemodialysis*. Hematol Oncol. [PMID 40684380](https://pubmed.ncbi.nlm.nih.gov/40684380/)
- Wen J et al. (2024). *Pharmacokinetics of Dexamethasone in Children and Adolescents with Obesity*. J Clin Pharmacol. [PMID 39120865](https://pubmed.ncbi.nlm.nih.gov/39120865/)
- Nakade S et al. (2008). *Population pharmacokinetics of aprepitant and dexamethasone in the prevention of chemotherapy-induced nausea and vomiting*. Cancer Chemother Pharmacol. [PMID 18317761](https://pubmed.ncbi.nlm.nih.gov/18317761/)
- Lennernäs H (2003). *Clinical pharmacokinetics of atorvastatin*. Clin Pharmacokinet. [PMID 14531725](https://pubmed.ncbi.nlm.nih.gov/14531725/)
- Schachter M (2005). *Chemical, pharmacokinetic and pharmacodynamic properties of statins: an update*. Fundam Clin Pharmacol. [PMID 15660968](https://pubmed.ncbi.nlm.nih.gov/15660968/)
- Tsamandouras N et al. (2017). *Modelling of atorvastatin pharmacokinetics and the identification of the effect of a BCRP polymorphism in the Japanese population*. Pharmacogenet Genomics. [PMID 27787353](https://pubmed.ncbi.nlm.nih.gov/27787353/)
- Liu Y et al. (2025). *Population Pharmacokinetics of Tranexamic Acid in Chinese Population Undergoing Cardiac Surgery with Cardiopulmonary Bypass*. Drug Des Devel Ther. [PMID 40453206](https://pubmed.ncbi.nlm.nih.gov/40453206/)
- Li S et al. (2021). *Population pharmacokinetics and pharmacodynamics of Tranexamic acid in women undergoing caesarean delivery*. Br J Clin Pharmacol. [PMID 33576009](https://pubmed.ncbi.nlm.nih.gov/33576009/)
- González Osuna A et al. (2022). *Population Pharmacokinetics of Intra-articular and Intravenous Administration of Tranexamic Acid in Patients Undergoing Total Knee Replacement*. Clin Pharmacokinet. [PMID 34255299](https://pubmed.ncbi.nlm.nih.gov/34255299/)
- Cirincione B et al. (2018). *Population Pharmacokinetics of Apixaban in Subjects With Nonvalvular Atrial Fibrillation*. CPT Pharmacometrics Syst Pharmacol. [PMID 30259707](https://pubmed.ncbi.nlm.nih.gov/30259707/)
- Byon W et al. (2017). *Population Pharmacokinetics, Pharmacodynamics, and Exploratory Exposure-Response Analyses of Apixaban in Subjects Treated for Venous Thromboembolism*. CPT Pharmacometrics Syst Pharmacol. [PMID 28547774](https://pubmed.ncbi.nlm.nih.gov/28547774/)
- Ueshima S et al. (2018). *Population pharmacokinetics and pharmacogenomics of apixaban in Japanese adult patients with atrial fibrillation*. Br J Clin Pharmacol. [PMID 29457840](https://pubmed.ncbi.nlm.nih.gov/29457840/)

## Outcomes & Prognosis

- Martinez-Perez R et al. (2022). *Third Ventricle Volume Predicts Functional Outcome in Chronic Subdural Hematoma*. Acta Neurol Scand. [PMID 34716574](https://pubmed.ncbi.nlm.nih.gov/34716574/)
- Knuutinen O et al. (2025). *Surgical Delay and Functional Outcome After Surgery for Chronic Subdural Hematoma*. World Neurosurg. [PMID 40024327](https://pubmed.ncbi.nlm.nih.gov/40024327/)
- Maldaner N et al. (2019). *Predicting Functional Impairment in patients with chronic subdural hematoma treated with burr hole Trepanation-The FIT-score*. Clin Neurol Neurosurg. [PMID 31125897](https://pubmed.ncbi.nlm.nih.gov/31125897/)
- Carlisi E et al. (2017). *Postoperative rehabilitation for chronic subdural hematoma in the elderly. An observational study focusing on balance, ambulation and discharge destination*. Eur J Phys Rehabil Med. [PMID 27145219](https://pubmed.ncbi.nlm.nih.gov/27145219/)
- Rauhala M et al. (2020). *Long-term excess mortality after chronic subdural hematoma*. Acta Neurochir (Wien). [PMID 32146525](https://pubmed.ncbi.nlm.nih.gov/32146525/)
- Petutschnigg T et al. (2026). *Long-Term Mortality, Cognition, and Quality of Life After Chronic Subdural Hematoma Surgery*. JAMA Neurol. [PMID 41973454](https://pubmed.ncbi.nlm.nih.gov/41973454/)
- Blaauw J et al. (2022). *Mortality after chronic subdural hematoma is associated with frailty*. Acta Neurochir (Wien). [PMID 36173514](https://pubmed.ncbi.nlm.nih.gov/36173514/)
- Zolin A et al. (2024). *Association of liver fibrosis with cognitive decline in Parkinson's disease*. J Clin Neurosci. [PMID 37976909](https://pubmed.ncbi.nlm.nih.gov/37976909/)
- Parikh NS et al. (2024). *Cognitive impairment and liver fibrosis in non-alcoholic fatty liver disease*. BMJ Neurol Open. [PMID 38268753](https://pubmed.ncbi.nlm.nih.gov/38268753/)
- Wahbeh F et al. (2024). *Impact of tobacco smoking on disease-specific outcomes in common neurological disorders: A scoping review*. J Clin Neurosci. [PMID 38428126](https://pubmed.ncbi.nlm.nih.gov/38428126/)

## QSP Methodology

- Scheuher B et al. (2024). *Towards a platform quantitative systems pharmacology (QSP) model for preclinical to clinical translation of antibody drug conjugates (ADCs)*. J Pharmacokinet Pharmacodyn. [PMID 37787918](https://pubmed.ncbi.nlm.nih.gov/37787918/)
- Gaweda AE et al. (2021). *Development of a quantitative systems pharmacology model of chronic kidney disease: metabolic bone disorder*. Am J Physiol Renal Physiol. [PMID 33308018](https://pubmed.ncbi.nlm.nih.gov/33308018/)
- Li X et al. (2023). *Combining network pharmacology, molecular docking, molecular dynamics simulation, and experimental verification to examine the efficacy and immunoregulation mechanism of FHB granules on vitiligo*. Front Immunol. [PMID 37575231](https://pubmed.ncbi.nlm.nih.gov/37575231/)
- Elmokadem A et al. (2019). *Quantitative Systems Pharmacology and Physiologically-Based Pharmacokinetic Modeling With mrgsolve: A Hands-On Tutorial*. CPT Pharmacometrics Syst Pharmacol. [PMID 31652028](https://pubmed.ncbi.nlm.nih.gov/31652028/)
- Lu T et al. (2024). *gPKPDviz: A flexible R shiny tool for pharmacokinetic/pharmacodynamic simulations using mrgsolve*. CPT Pharmacometrics Syst Pharmacol. [PMID 38082557](https://pubmed.ncbi.nlm.nih.gov/38082557/)
- Choules MP et al. (2024). *Physiologically based pharmacokinetic model to predict drug-drug interactions with the antibody-drug conjugate enfortumab vedotin*. J Pharmacokinet Pharmacodyn. [PMID 37632598](https://pubmed.ncbi.nlm.nih.gov/37632598/)
- Li Q et al. (2024). *Physiologically-Based Pharmacokinetic Modeling and Dosing Optimization of Cefotaxime in Preterm and Term Neonates*. J Pharm Sci. [PMID 38460573](https://pubmed.ncbi.nlm.nih.gov/38460573/)

---

## What the Literature Does Not Supply

The values below were not measured directly in the literature but were **estimated by this model to fit observed clinical outcomes**. They must be mentioned alongside any citation of the model.

- Neomembrane hydraulic conductivity `LP_IMM`/`LP_MAT` and reflection coefficients `SIG_IMM`/`SIG_MAT` — the qualitative difference between immature and mature vessels is well established, but no Starling coefficients have been directly measured in cSDH neomembranes.
- Clearance capacity via dural lymphatics `K_LYMPH` — the existence of the pathway is established, but its quantitative capacity is unknown.
- The Hill exponent of the haem switch `HB_HILL = 2` and threshold `HB50` — the key assumption that creates bistability, and the structural claim in this model most in need of validation.
- The atrophy-dependence of brain re-expansion rate `K_REEXP0·exp(−V_RES/LAM_ATR)` — a clinically well-known phenomenon, but the functional form is this model's choice.
- The cavity-patency term `f_pat` (D50_PAT, N_PAT) — governs how fast the drained space refills, but there is no direct empirical evidence for it.
