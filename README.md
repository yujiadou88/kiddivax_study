# Hong Kong Influenza Vaccination study

## 1. Study overview

While vaccination is known to be effective in reducing influenza-related morbidity in school-age children it is not currently recommended in Hong Kong because few infections in this age group above 2 years result in significant morbidity or mortality. Also, the supply of vaccines is limited and indirect costs of vaccination are high.

However, evidence from randomized trials, observational studies and modeling studies has accumulated which suggests that the vaccination of school-age children can provide substantial indirect benefits to the rest of the community, particularly to older adults who are most at risk for the complications associated with influenza. This approach to vaccination is known as the "transmission-limiting" strategy, in contrast to the current "morbidity-limiting" approach where those individuals at highest risk of complications are targeted to receive the vaccine.

While the indirect benefits of vaccinating school-age children against influenza have been extensively studied at the community level, relatively few studies have investigated whether vaccination of school-age children can lead to indirect benefits at the household level, to the parents or siblings of vaccinees.

In 2008-09 (pilot) and 2009-10 (main study) we conducted a double-blind randomized controlled trial of households in Hong Kong, where one child aged 6-17 from each household was randomized to receive either one dose of trivalent inactivated seasonal influenza vaccine, or saline placebo. Each household was followed up for around a year, with serum drawn from each household member before and after the study, and also from some participants during the study. Subjects and household contacts were asked to keep daily symptom diaries, and during episodes of acute respiratory illness in any household member we arranged to collect nose and throat swabs from all household members regardless of illness for laboratory confirmation of influenza infection.

## 2. Raw data

The latest version of the study year 1 (pilot) data (09-2008 to 12-2009) are available to download as a zip file here:

-   [KiddivaxPilot.zip](data/KiddivaxPilot.zip)

This version of the dataset covers the incidence of infections and illnesses in household members, and serologic and virologic laboratory results.

The latest version of the study year 2 (main) data (09-2009 to 12-2010) are available to download as a zip file here:

-   [KiddivaxMain.zip](data/KiddivaxMain.zip)

We provide our data under the [Open Data Commons Public Domain Dedication and License](http://www.opendatacommons.org/odc-public-domain-dedication-and-licence/), which is a version of open access for data. Under this licence we reserve no rights: there are no restrictions on use of our data, and no requirement to cite our work or this website. However we would anticipate that for academic purposes the standard practice of referencing sources would apply. We would like to hear from researchers who are using our data and we would be keen to work together on analyses.

## 3. Main findings

The primary findings from the 2008-09 pilot study were published by [Cowling et al. (2010, Clin Infect Dis)](http://dx.doi.org/10.1086/657311). In brief, the pilot study achieved its aim of confirming the feasibility of the main study design, provided some information about the incidence of influenza infections and illnesses, and provided primary evidence of temporary immunity. We found no evidence of indirect benefits although with only 119 households our pilot study was underpowered to detect even fairly large indirect effects.

Results presented in that manuscript are reproduced in the following scripts which can be run in the free (and [increasingly popular](http://www.nytimes.com/2009/01/07/technology/business-computing/07program.html)) statistical software package [R](http://www.r-project.org):

-   [scripts](KiddivaxPilot_scripts/dataframe.r) to reformat the raw data (used in some of the other scripts here) and [scripts](kiddivaxPilot_scripts/dataframe_cross.r) to reformat the raw data and adjust for cross reactivity.

-   [Table 1](KiddivaxPilot_scripts/Table_1.r)

-   [Table 2](KiddivaxPilot_scripts/Table_2.r)

-   [Table 3](KiddivaxPilot_scripts/Table_3.r)

-   [Table 4](KiddivaxPilot_scripts/Table_4.r)

-   [Table 5](KiddivaxPilot_scripts/Table_5.r)

-   [Table 6](KiddivaxPilot_scripts/Table_6.r)

-   [Table 7](KiddivaxPilot_scripts/Table_7.r)

-   [Table 8](KiddivaxPilot_scripts/Table_8.r)

-   [Figure 1](KiddivaxPilot_scripts/Figure_1.r)

-   [Figure S1](KiddivaxPilot_scripts/Figure_S1.r)

-   [Figure S2](KiddivaxPilot_scripts/Figure_S2.r)

-   [Figure S3](KiddivaxPilot_scripts/Figure_S3.r)

-   [Figure S4](KiddivaxPilot_scripts/Figure_S4.r)

## 4. Transmissibility of influenza in households

We analyzed the final outbreak size distributions of pH1N1, sH1N1, sH3N2 infections identified in paired sera collected from 117 Hong Kong households in April and in August-October 2009. The results were published by [Klick et al. (2011, Epidemiology)](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3206962/). In conclusion, pandemic and seasonal influenza A viruses had similar age-specific transmissibility in a cohort of initially uninfected households, after adjustment for baseline immunity.

Results described in Klick et al. (2011) are reproduced in the following scripts which can be run in [R](http://www.r-project.org):

-   [functions](kiddivaxPilot_Transmission_scripts/MCMC_function.r) of the MCMC procedure (used in some of the other scripts here).

-   [Table 1](KiddivaxPilot_Transmission_scripts/Table_1.r); [Table 2](KiddivaxPilot_Transmission_scripts/Table_2.r); [eTable 2](KiddivaxPilot_Transmission_scripts/eTable_2.r); [eTable 3](KiddivaxPilot_Transmission_scripts/eTable_3.r); [eTable 4](KiddivaxPilot_Transmission_scripts/eTable_4.r); [eTable 5](KiddivaxPilot_Transmission_scripts/eTable_5.r).

-   [Figure 1-3](KiddivaxPilot_Transmission_scripts/Figure_1to3.r).

## 5. Increased risk of noninfluenza respiratory virus infections

We randomized 115 children to TIV or placebo group. TIV recipients had an increased risk of noninfluenza respiratory virus infections over the following nine months. While being protected against influenza, TIV recipients may lack of temporary non-specific immunity that protected against other respiratory viruses. The results were published by [Cowling et al. (2012, Clin Infect Dis)](http://cid.oxfordjournals.org/content/early/2012/03/13/cid.cis307.short).

Results described in the manuscript are reproduced in the following scripts which can be run in [R](http://www.r-project.org):

-   [Table 1](KiddivaxPilot_CID_respiratory_scripts/Table_1.r); [Table 2](KiddivaxPilot_CID_respiratory_scripts/Table_2.r); [Table 3](KiddivaxPilot_CID_respiratory_scripts/Table_3.r); [Appendix Table 1](KiddivaxPilot_CID_respiratory_scripts/Appendix_Table_1.r).

-   [Figure 1](KiddivaxPilot_CID_respiratory_scripts/Figure_1.r); [Appendix Figure 1](KiddivaxPilot_CID_respiratory_scripts/Appendix_Figure_1.r).

## 6. Main findings

The primary findings from the 2009-10 main study were published by [Cowling et al. (2012, Clin Infect Dis)](http://cid.oxfordjournals.org/content/55/5/695). In the main study, we recruited 796 children as well as their household members. Pandemic A(H1N1) circulated at the time of vaccination and for a short time afterward with no substantial seasonal influenza activity during that period. We found that seasonal TIV prevented pandemic influenza A(H1N1) and influenza B infections in children.

Results described in the manuscript are reproduced in the following scripts which can be run in [R](http://www.r-project.org):

-   [scripts](KiddivaxMain_scripts/dataframe.r) to reformat the raw data (used in some of the other scripts here).

-   [Table 1](KiddivaxMain_scripts/Table_1.r); [Table 2](KiddivaxMain_scripts/Table_2.r); [Table 3](KiddivaxMain_scripts/Table_3.r); [Table S1](KiddivaxMain_scripts/Table_S1.r); [Table S2](KiddivaxMain_scripts/Table_S2.r); [Table S3](KiddivaxMain_scripts/Table_S3.r); [Table S4](KiddivaxMain_scripts/Table_S4.r).

-   [Figure 2](KiddivaxMain_scripts/Figure_2.r).

## 7. Humoral antibody response after receipt of TIV

Among 119 children recruited in the pilot study, 64 children rejoined the main study and allocation of TIV/placebo in year 2 was independent of allocation in year 1. [Ng et al. (2012, Pediatr Infect Dis J)](http://journals.lww.com/pidj/Fulltext/2012/09000/Humoral_Antibody_Response_After_Receipt_of.24.aspx) studied those 64 subjects and suggested that humoral antibody response to TIV may be lower in children receiving repeated vaccination, but receipt of TIV induced seroprotection in most subjects.

Results described in the manuscript are reproduced in the following scripts which can be run in [R](http://www.r-project.org):

-   [scripts](KiddivaxMain_PIDJ_antibody_scripts/dataframe.r) to reformat the raw data (used in some of the other scripts here).

-   [Table 1](KiddivaxMain_PIDJ_antibody_scripts/Table_1.r); [Table 2](KiddivaxMain_PIDJ_antibody_scripts/Table_2.r); [Table 3](KiddivaxMain_PIDJ_antibody_scripts/Table_3.r).

-   [Figure 1](KiddivaxMain_PIDJ_antibody_scripts/Figure_1.r).

## 8. The effect of age and recent vaccination history on the efficacy of 2009-10 seasonal TIV

There is some evidence that annual vaccination of trivalent inactivated influenza vaccine (TIV) may lead to reduced vaccine immunogenicity but evidence is lacking on whether vaccine efficacy is affected by prior vaccination history. [Ng et al. (2013, PLoS ONE)](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0059077) investigated on the 796 children recruited in the main study, whose influenza vaccination history in the two preceding years was recorded. Ng concluded that prior vaccination was associated with lower antibody responses to TIV against seasonal influenza A vaccine viruses, but higher responses to the same lineage of influenza B virus.

Results described in the manuscript are reproduced in the following scripts which can be run in [R](http://www.r-project.org):

-   [scrits](KiddivaxMain_PLoSONE_scripts/dataframe.r) to reformat the raw data (used in some of the other scripts here).

-   [Table 2](KiddivaxMain_PLoSONE_scripts/Table_2.r); [Table 3](KiddivaxMain_PLoSONE_scripts/Table_3.r); [Table S1](KiddivaxMain_PLoSONE_scripts/Table_S1.r); [Table S2](KiddivaxMain_PLoSONE_scripts/Table_S2.r); [Table S3](KiddivaxMain_PLoSONE_scripts/Table_S3.r); [Table S4](KiddivaxMain_PLoSONE_scripts/Table_S4.r).

-   [Figure 1](KiddivaxMain_PLoSONE_scripts/Figure_1.r); [Figure S2](KiddivaxMain_PLoSONE_scripts/Figure_S2.r).

## 9. The association between antibody titers and protection against influenza

Antibody titers measured by hemagglutination inhibition (HAI) correlate with protection against influenza virus infection and are used to specify criteria for vaccine licensure. [Ng et al. (2013, JID)](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3778972/) investigated on the 773 children recruited in the main study with HAI titers available, and found that HAI titers of 1:40 against A(H1N1)pdm09 and B(Victoria lineage) were associated with 48% and 55% protection against PCR-confirmed infection with each strain.

Results described in the manuscript are reproduced in the following scripts which can be run in [R](http://www.r-project.org):

-   [source 1](KiddivaxMain_JID_titerprot_scripts/source_1.r); [source 2](KiddivaxMain_JID_titerprot_scripts/source_2.r); [source 3](KiddivaxMain_JID_titerprot_scripts/source_3.r); [source 4](KiddivaxMain_JID_titerprot_scripts/source_4.r); [source 5](KiddivaxMain_JID_titerprot_scripts/source_5.r); [source 6](KiddivaxMain_JID_titerprot_scripts/source_6.r) (used in some of the other scripts here).

-   [Figure 1](KiddivaxMain_JID_titerprot_scripts/Figure_1.r); [Figure 2](KiddivaxMain_JID_titerprot_scripts/Figure_2.r); [Appendix Figure 1](KiddivaxMain_JID_titerprot_scripts/Appendix_Figure_1.r);

-   [Appendix Table 1](KiddivaxMain_JID_titerprot_scripts/Appendix_Table_1.r); [Appendix Table 2](KiddivaxMain_JID_titerprot_scripts/Appendix_Table_2.r).

## 10. Further work

We have a series of further analyses underway using the data from our study. More details will follow later.

## Authors and investigator

The principal investigator of this study is [Ben Cowling](https://sph.hku.hk/en/Biography/Cowling-Benjamin-John). The data were uploaded by [Ben Cowling](https://sph.hku.hk/en/Biography/Cowling-Benjamin-John), and the scripts were written by [Ben Cowling](https://sph.hku.hk/en/Biography/Cowling-Benjamin-John), Vicky Fang, Brendan Klick, and Sophia Ng.

## A comment on reproducible research

We fully support the [increasing calls](http://dx.doi.org/10.1097/EDE.0b013e318196784a) from the academic community for [epidemiologic analyses to be reproducible](http://dx.doi.org/10.1093/aje/kwj093 "Peng et al., 2006, AJE"), and [raw data from randomized controlled trials to be published](http://dx.doi.org/10.1186/1745-6215-10-17 "Hrynaszkiewicz and Altman, 2009, Trials"), as a part of the wider scientific endeavour to replicate results. Another example of this recommendation is in the [Good Practice Guide for Quantitative Veterinary Epidemiology](http://www.qve-goodpracticeguide.org.uk/guide#TOC-2.4-Inputs). Here we have published the raw *anonymised* data from our HK NPI pilot study, and will later release the data from our main study. We have also published scripts which allow the analyses in our published papers to be reproduced.

Thousands of local people have given their time, and their families', as part of their participation in our studies, all in the expectation that our research studies will add to medical and scientific knowledge. Participants should also expect that we will make the best possible use of the information that we have collected about them. It would be difficult to argue that facilitating best use of the data by the research community need not involve releasing raw data.

Publication of anonymised raw data has been approved by our local IRB and funding sources, and participants were advised that anonymised data would be published during the informed consent process. We anticipate that release of the raw data will:

-   Promote reproducibility of our results.

-   Allow other investigators to conduct their own analyses on our data.

-   Allow other investigators to compare our data with theirs, for example to explore similarities and differences between research findings.

-   Allow other investigators to prepare and plan for their own studies.

## Publications

1.  Cowling BJ, Ng S, Ma ESK, Cheng CKY, Wai W, Fang VJ, Chan KH, Ip DKM, Chiu SS, Peiris JSM, and Leung GM. Protective efficacy of seasonal influenza vaccination against seasonal and pandemic influenza virus infection during 2009 in Hong Kong. *Clinical Infectious Diseases*, 2010; 51(12):1370-9. [[link]](http://dx.doi.org/10.1086/657311) [[PubMed]](http://www.ncbi.nlm.nih.gov/pubmed/21067351).

2.  Klick B, Nishiura H, Ng S, Fang VJ, Leung GM, Peiris JSM, and Cowling BJ. Transmissibility of seasonal and pandemic influenza in a cohort of households in Hong Kong in 2009. *Epidemiology*, 2011; 22(6):793-6. [[link]](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3206962) [[PubMed]](http://www.ncbi.nlm.nih.gov/pubmed/21878814).

3.  Cowling BJ, Fang VJ, Nishiura H, Chan KH, Ng S, Ip DKM, Chiu SS, Leung GM and Peiris JSM. Increased risk of noninfluenza respiratory virus infections associated with receipt of inactivated influenza vaccine. *Clinical Infectious Diseases*, 2012; 54(12):1778-83. [[link]](http://cid.oxfordjournals.org/content/54/12/1778) [[PubMed]](http://www.ncbi.nlm.nih.gov/pubmed/22423139).

4.  Cowling BJ, Ng S, Ma ES, Fang VJ, So HC, Wai W, Cheng KY, Wong JY, Chan KH, Ip DKM, Chiu SS, Peiris JSM, and Leung GM. Protective Efficacy Against Pandemic Influenza of Seasonal Influenza Vaccination in Children in Hong Kong: A Randomized Controlled Trial. *Clinical Infectious Diseases*, 2012; 55(5):695-702. [[link]](http://cid.oxfordjournals.org/content/55/5/695) [[PubMed]](http://www.ncbi.nlm.nih.gov/pubmed/22670050).

5.  Ng S, Fang VJ, Ip DKM, Chiu SS, Leung GM, Peiris JSM, and Cowling BJ. Humoral Antibody Response after Receipt of Inactivated Seasonal Influenza Vaccinations One Year Apart in Children. *The Pediatric Infectious Disease Journal*, 2012; 31(9):964-9. [[link]](http://journals.lww.com/pidj/Fulltext/2012/09000/Humoral_Antibody_Response_After_Receipt_of.24.aspx) [[PubMed]](http://www.ncbi.nlm.nih.gov/pubmed/22683675).

6.  Ng S, Ip DKM, Fang VJ, Chan KH, Chiu SS, Leung GM, Peiris JSM, and Cowling BJ. The effect of age and recent influenza vaccination history on the immunogenicity and efficacy of 2009-10 seasonal trivalent inactivated influenza vaccination in children. *PLoS ONE*, 2013; 8(3):e59077. [[link]](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0059077). [[PubMed]](http://www.ncbi.nlm.nih.gov/pubmed/23554974).

7.  Ng S, Fang VJ, Ip DKM, Chan KH, Leung GM, Peiris JSM, and Cowling BJ. Estimation of the association between antibody titers and protection against confirmed influenza virus infection in children. *The Journal of Infectious Diseases*, 2013; 208(8):1320-4. [[link]](http://jid.oxfordjournals.org/content/208/8/1320.long) [[PubMed]](http://www.ncbi.nlm.nih.gov/pubmed/23908481).

## Acknowledgements

The pilot study (2008-09) was funded by the Area of Excellence Scheme of the Hong Kong University Grants Committee (grant no. AoE/M-12/06) and the Research Fund for the Control of Infectious Disease, Food and Health Bureau, Government of the Hong Kong SAR (grant no. PHE-02).

The main study (2009-10) was funded by the Research Fund for the Control of Infectious Disease, Food and Health Bureau, Government of the Hong Kong SAR (grant no. CHP-CE-03).

Further analyses are supported by the Area of Excellence Scheme of the Hong Kong University Grants Committee (grant no. AoE/M-12/06), and the Harvard Center for Communicable Disease Dynamcs from the US National Institutes of Health Models of Infectious Disease Agent Study program (grant no. 1 U54 GM088558).

------------------------------------------------------------------------

[![Creative Commons License](https://i.creativecommons.org/l/by/3.0/80x15.png)](http://creativecommons.org/licenses/by/3.0/) This work is licensed under a [Creative Commons Attribution 3.0 Unported License](http://creativecommons.org/licenses/by/3.0/). It was written by [Ben Cowling](http://www.hku.hk/cmd/staff/bio/cowling.htm) (email).
