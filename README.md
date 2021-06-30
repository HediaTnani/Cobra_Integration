
# Identification of novel drug targets in Leishmania major using in-silico metabolic coupling

Cutaneous Leishmaniasis is a debilitating disease for approximately 12 million
people living mainly within the tropics. It is caused by different species of the protozoan
Leishmania parasites among which Leishmania major. The complex interactions
between Leishmania and macrophages is central to the outcome of the disease (Pimenta
et al. 1997). The lack of an effective vaccine and vector control program makes the use
of chemotherapy the primary tool against leishmaniasis. However, the emergence of
resistance against available chemotherapies is a major concern and make the discovery
of new drug targets and the development of less toxic drugs compulsory.

# Methodology

## 1- Macrophage model
We built a non-infected mouse macrophage model using GIMME algorithm
(Becker and Palsson 2008) from the COBRA Toolbox v3.0 (Heirendt et al. 2019) by
using the gene expression data of non-infected mice (Rabhi et al. 2012), and the latest
mouse metabolic model iMM1865 (Khodaee et al. 2020).


1. Download the series matrix format of gene expression data from GEO Omnibus

2. The following R code was used for data preprocessing:

Code: Data_prep_gimme.R 

4. The expression dataset with IDs corresponding to iMM1865 were imported into MATLAB having COBRA Toolbox v3.0. The following code was run to map the expression
to the model using GIMME algorithm. 

Code: gimme.m
