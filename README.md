# Harmonize

Created by: Seehyun Park

Creation date: 01 Jan 2023  

https://cesp.inserm.fr/en/equipe/exposome-and-heredity  


## A tutorial for performing harmonization of GWAS summary statistics 

## Contents
- [Note](#notes)
- [Purpose](#purpose)
- [Steps](#steps)

        
## <a id="notes" /> Notes
The use of external Summary Statistics in genome-wide association study (GWAS) can significantly increase the size and diversity of the sample, increasing the power to detect association analysis. However, due to batch effects, genotyping errors, and the use of different genotyping platforms, the aggregation of multiple GWAS summary statistics can be quite challenging and difficult. If these GWAS summary statistics are not carefully quality controlled, the incorrect results might be derived when performing meta-analysis. 


## <a id="purpose" /> Purpose
The purpose of this package is to provide the data harmonization pipeline that allows users to check the quality of GWAS summary statistics before performing a meta-analysis. 

## <a id="steps" /> Steps
**Standardize** 
The standardize module scans throught the GWAS summary statistics and remove the SNPs according to the following criteria
- Null: Null or NA values are genotyped on one of the beta, effect allele frequency or p-value
- Weird: Genotyped are not based on the combination of 'ATCG'
- Duplicate: The same SNP is expressed in duplication
- Palindromic: The palindromic SNP, such that the alleles on the forward strand are the same as on the reverse strand (A/T on forward is T/A on the reverse). However 






![Image4](https://user-images.githubusercontent.com/60399683/208106098-36278287-dc55-47a4-adc5-5542a5ca00f6.png)

In fact, these SNPs associated with the exposure are used as proxy to determine the exposure. This method is analogous to randomized controlled trials because in MR, the genetic alleles of these SNPs are randomly distributed at birth (in a rct, it's the treatment randomly assigned).  

In a two sample MR, the association *(i)* between the SNPs and the exposure and *(ii)* between the SNPs and the outcome (needed to estimate the causal estimate between the exposure and the outcome) come from two different and **INDEPENDANT** samples (without overlapping !!) 

## <a id="dataprepa" /> Data preparation

### <a id="expo" /> Exposures file

In this study, we assess three exposures: smoking, alcohol and coffee consumption. The IVs (SNPs associated to these exposures, gwas significant 5e-8) come from the largest GWAS available in the litterature at the time of the analysis (2021). (Liu 2019 for smoking and alcohol, Zhong 2019 for coffee).
Report these association in a .txt file (let call this file "exposures_noclump.txt") with a header like this : 
|Phenotype|SNP|CHR|POS|EA|BA|EAF|BETA|SE|PVAL|N|
|---------|---|---|---|--|--|---|----|--|----|-|

One of the column (the first for exemple) must indicate the exposure of interest for the selected SNPs and must be called "Phenotype" so that the R package understands that it concerns the different exposures studied. (EA means Effect Allele, BA means Base Allele, i.e. non effect allele, N: number total of sample included in the analysis)

#### <a id="clump" /> Clump the data 

First of all, all the SNPs have to be independant. The recommand threshold for linkage desequilibrium (LD) is very low in MR: r²=0.001. 
Lets do it: 
```
library(TwoSampleMR)

exposure = read.table("exposures_noclump.txt", header=T, sep="\t")

expo = read_exposure_data(
  filename = "exposures.txt",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "EA",
  other_allele_col = "BA",
  eaf_col = "EAF",
  pval_col = "PVAL")
  
expo2 = clump_data(expo, pop="EUR")

expo_clump = merge(exposure, expo2[,c(1,2)], by="SNP")
expo_clump = expo_clump[order(expo_clump$Pheno, expo_clump$CHR, expo_clump$POS), c(2,1,3:11)]

```

#### <a id="beta" /> Make the beta of the exposure positive

One method used in our analysis is the "MR-Egger" method. When use this method, the betas of the association between the SNP and the exposure have to be positive (maybe in the current version of the package it's not necessary, but here we do this step to be sure). You can do this using this script:

```
expo_clump = read.table("0_data/exposures_0.001.txt", header=T, sep="\t")

expo_clump$EAb = ifelse(expo_clump$BETA > 0, expo_clump$EA, expo_clump$BA)
expo_clump$BAb = ifelse(expo_clump$BETA > 0, expo_clump$BA, expo_clump$EA)
expo_clump$EAFb = ifelse(expo_clump$BETA > 0, expo_clump$EAF, 1-expo_clump$EAF)
expo_clump$BETAb = abs(expo_clump$BETA)
```

#### <a id="prox" /> Looking for proxies?

**TIPS**: At this point, you can check if all the SNPs used as IVs are available in your GWAS on the outcome of interest. In this study, we studied the association with Parkinson's disease (PD) and we used a meta-analysis on European individuals from the inernational Courage-PD consortium. 
You can search proxies in high LD (for exemple r²>0.8) in [LDLink](https://ldlink.nci.nih.gov/?tab=ldproxy) or [SNiPA](https://snipa.helmholtz-muenchen.de/snipa3/). We advice to use the SNP with the higher LD to have a very good correspondance of the alleles (the [LDpair](https://ldlink.nci.nih.gov/?tab=ldpair) tool on LDLink is very useful to find the correspondance of the alleles).

Here, we find that the SNP rs6050446 used as IV for the "Smoking" exposure was not present in our Courage-PD meta-analysis, but one of its proxy, the rs117495226
in very high LD with it (r²=0.9, D'=1) was present in Courage-PD. So we used it instead of the initial IV. 
**BE CAREFUL** : If the SNP is palindromic, you must check the allele frequency to determine the correspondance of the alleles. If the MAF is too close to 0.50, you can't determine the correspondance, so you can't use the SNP. 

Here the script to incorporate this new information (you also can do it manually using Excel):

```
expo_clump$Proxy_SNP=as.character(NA)
expo_clump$Proxy_EA=as.character(NA)
expo_clump$Proxy_BA=as.character(NA)
expo_clump$Proxy_EAF=as.numeric(NA)
expo_clump$Proxy_DPRIME=as.numeric(NA)
expo_clump$Proxy_R2=as.numeric(NA)
expo_clump$Proxy_POS=as.numeric(NA)

expo_clump$Proxy_SNP = ifelse(expo_clump$SNP=="rs6050446", "rs117495226", expo_clump$Proxy_SNP)
expo_clump$Proxy_EA = ifelse(expo_clump$Proxy_SNP=="rs117495226","G", expo_clump$Proxy_EA)
expo_clump$Proxy_BA = ifelse(expo_clump$Proxy_SNP=="rs117495226","C", expo_clump$Proxy_BA)
expo_clump$Proxy_EAF = ifelse(expo_clump$Proxy_SNP=="rs117495226",0.95, expo_clump$Proxy_EAF)
expo_clump$Proxy_DPRIME = ifelse(expo_clump$Proxy_SNP=="rs117495226",1, expo_clump$Proxy_DPRIME)
expo_clump$Proxy_R2 = ifelse(expo_clump$Proxy_SNP=="rs117495226",0.9, expo_clump$Proxy_R2)
expo_clump$Proxy_POS = ifelse(expo_clump$Proxy_SNP=="rs117495226",25565515, expo_clump$Proxy_POS)

# Final variable after taking into account the proxies (variable_f for variable_final)

expo_clump$SNP_f = ifelse(is.na(expo_clump$Proxy_SNP), expo_clump$SNP, expo_clump$Proxy_SNP)
expo_clump$CHR_f = expo_clump$CHR
expo_clump$POS_f = ifelse(is.na(expo_clump$Proxy_SNP), expo_clump$POS, expo_clump$Proxy_POS)
expo_clump$EA_f = ifelse(is.na(expo_clump$Proxy_SNP), expo_clump$EAb, expo_clump$Proxy_EA)
expo_clump$BA_f = ifelse(is.na(expo_clump$Proxy_SNP), expo_clump$BAb, expo_clump$Proxy_BA)
expo_clump$EAF_f = ifelse(is.na(expo_clump$Proxy_SNP), expo_clump$EAFb, expo_clump$Proxy_EAF)
expo_clump$BETA_f = expo_clump$BETAb
expo_clump$SE_f = expo_clump$SE
expo_clump$PVAL_f = expo_clump$PVAL

write.table(expo_clump, "exposures_0.001.txt", quote=F, row.names=F, sep="\t")
```
You can find the "exposure_0.001.txt" file in the 0_data folder here.

So, following this script, the final variables used for the MR are: SNP_f, CHR_f, POS_f, EA_f, BA_f, EAF_f, BETA_f, SE_f, PVAL_f

### <a id="outcome" /> Outcome file

As with the exposures, you will need to create a file listing the associations between the snps used as VIs (definied as SNP_f in the exposures_0.001.txt file) in your analysis and your outcome of interest. Here, these associations come from the Courage-PD meta-analysis as mentioned above. 
You can extract these SNPs using their coordinates (!!! Be carreful of the build used - hg19 or hg38 for exemple).  

For exemple, 
1. In the whole GWAS (gwas.txt), retain the number of column of the chromosome and position (with this exemple, chr is column one and position the column two):
   |CHR|POS|SNP|EA|BA|BETA|SE|PVAL|I2|Q.VALUE|N_STUDIES|N|EAF|
   |---|---|---|--|--|----|--|----|--|-------|---------|-|---|

2. Create a list_snp.txt file containing the list of the SNPs used as IVs you want to extract from the whole gwas, like this:
   |CHR|POS|SNP|
   |---|---|---|

3. You can use this AWK command to extract them:

```
# If the names of the "chromosome" and "position" column in both files are the same (like here), use this command:
awk 'NR==FNR {a[$1":"$2]=$1; next} $1":"$2 in a {print}' list_snp.txt gwas.txt > extraction_snp.txt

# If the names are different (for exemple "CHR" and "POS" in the gwas.txt file and "CHROMOSOME" and "POSITION" for the list_snp.txt file, use this command:
head -1 gwas.txt > extraction_snp.txt
awk 'NR==FNR {a[$1":"$2]=$1; next} $1":"$2 in a {print}' list_snp.txt gwas.txt >> extraction_snp.txt

```

It is possible to assess several outcomes. For this, you can add the "Phenotype" column containing the name of the different outcomes (as it works for the exposures.txt file).
Here, the outcome file (available in the 0_data, called outcome.txt) looks like this:
|Phenotype|SNP|CHR|POS|EA|BA|BETA|SE|PVAL|I2|Q.VALUE|N_STUDIES|N|EAF|
|---------|---|---|---|--|--|----|--|----|--|-------|---------|-|---|

(Please note that in this study, you used only SNPs that were available in 17/23 Courage-PD studies (75%))

## <a id="2sampleMR" /> Performing Two Sample MR

### <a id="fun" /> Functions definition

In the 00_scripts repository, the 0_functions.R contains all the function automatically implemented to run the MR analyses. To change the paths, you have to change them in this file. 
Please run these functions in the first time to load them in the R environment. They will be explained in the Readme file. 

(You can use the commandArgs function to use arguments with Rscript)

### <a id="mainMR" /> Main analyses

After running the 00_scripts/0_functions.R, the command to run the MR analyses are described here, and also written in the 00_scripts/1_runMR.R script

In the main analysis, we used several methods to estimate the causal effect.
- The IVW, which is the principal analysis for MR if IVs>1 (if there is only one IV, you can use only the wald ratio method = beta(SNP-outcome)/beta(SNP-expo)
  The IVX fix and random are shown in the table of results, but it is recommand to use the random everytime (please see [this article](https://wellcomeopenresearch.org/articles/4-186#:~:text=The%20guidelines%20are%20divided%20into,%2C%20data%20presentation%2C%20and%20interpretation.), a very important article to understand the MR).
- The MR-Egger method (for horizontal directional pleiotropy)
- The Weighted median and mode based
- The MR-PRESSO
The main function to automatically write the table of results is called "MR_fun".

> MR_fun(name_exposure_file, name_outcome_file, strata) {}

- **name_exposure_file** is the name of the file with the associations between SNPs and exposures after clump (here: "exposures_0.001")
- **name_outcome_file** is the name of the file with the associations between SNPs and outcome (here: "outcome")
- **strata** is the name of the strata of the analysis (here: "All indiv")
(these files are in the 0_data folder)

```
list_presso = list("Smoking", "Alcohol", "Coffee")

MR_fun("exposures_0.001", "outcome", "All indiv")
```

Before running the MR_fun, you have to define the list of your exposures in the list_presso object to use them in the MR_presso function. Fill in only the exhibitions for which there are several VI (about 5 minimum). 

Running this function automatically create the Results_MR_exposures_0.001.xlsx file in the 1_results folder.
This file contains several columns:

|exposure|method|nsnp|b|se|OR|pval|pval_intercept|Q_pval|P_Glob|P_Dist|Outliers|beta_presso|se_presso|OR_presso|p_presso|
|--------|------|----|-|--|--|----|--------------|------|------|------|--------|-----------|---------|---------|--------|

Here the definition of these different columns:
- **exposure**: contains the name of the exposure assessed (here: Smoking, Alcohol and Coffee)
- **method**: contains the method of MR used (IVW fix and random, MR-Egger, Simple mode, Weighted mode, Weighted median). Please note that the MR presso results are not in this column (see below).
- **nsnp**: contains the number of SNPs used for the MR after removing SNPs not included in the outcome gwas, palindromic etc. (automaticaly done by the package)
- **b**: contains the beta of the MR association (*i.e* between the exposure and the outcome) estimated by the method definied in the method column
- **se**: same for standard error
- **OR**: same for exp(beta) and 95% confidence interval exp(beta)
- **pval**: same for pvalue
- **pvalue_intercept**: is the pvalue of the MR-Egger intercept. If the pvalue is significant (<0.05), it means that there is directional pleiotropy and that you have to take into account the estimate provided by the MR-Egger method, the one estimated by the IVW method is then biaised and not true. Otherwise, if the intercept of MR-Egger is not significant, it is not necessary to take into account the estimation made by the MR-Egger method.
- **Q_pval**: is the pvalue of the Cochran Q test for heterogeneity. If Q_pval<0.05, SNPs in the instrument are heterogenous. 
- The others columns are the results of the MR-PRESSO method. This is a 3 in 1 method:
  - **P_Glob**: first of all, the MR-Presso method estimates the global horizontal pleiotropy (*i.e* the global heterogeneity) in the instrument. This is the "Global test". If P_Glob<0.05, there is pleiotropy according to MR-Presso method. If MR_presso finds heterogeneity, it will try to identify outliers SNPs, and in this case, remove them. It is possible that MR_presso showed heterogenity (P_Glob<.05) but does not identfied outliers.
  - **Outliers**: the names of the outliers SNPs identified and removed are mentionned in this column
  - **beta_presso**, se_presso, OR_presso and p_presso are the estimation corrected by MR_Presso after removing these outliers
  - **P_Dist**: Finally, MR_Presso estimates if there is a significant difference between the original beta (before removing outliers) and the corrected beta (after removing outliers) by doing a "distortion test". If P_Dist<0.05, there is a signficant differance between both estimates.

For exemple, here are the results for the Alcohol and Smoking exposures:

| exposure | method          | nsnp | OR (95% CI)    | pval | pval_intercept | Q_pval | P_Glob | P_Dist | Outliers | OR_presso (95% CI)| p_presso |
| -------- | --------------- |------|----------------|------|----------------|--------|--------|--------|----------|-------------------|----------|
| Alcohol  | IVW fix         |      |0.72 (0.45-1.16)| 0.18 |                |        |        |        |          |                   |          |
| Alcohol  | IVW random      | 62   |0.68 (0.39-1.18)| 0.17 |                |  0.06  |  0.047 |  0.71  |rs9607814 |0.77 (0.46-1.29)   |   0.33   |
| Alcohol  | MR Egger        | 62   |0.70 (0.30-1.61)| 0.40 | 0.94           |  0.05  |        |        |          |                   |          | 
| Alcohol  | Simple mode     | 62   |0.32 (0.04-2.58)| 0.29 |                |        |        |        |          |                   |          |
| Alcohol  | Weighted median | 62   |0.86 (0.42-1.75)| 0.67 |                |        |        |        |          |                   |          |
| Alcohol  | Weighted mode   | 62   |0.86 (0.44-1.69)| 0.66 |                |        |        |        |          |                   |          |
| Smoking  | IVW fix         |      |0.79 (0.63-0.97)| 0.03 |                |        |        |        |          |                   |          |
| Smoking  | IVW random      | 182  |0.74 (0.60-0.93)| 0.01 |                |  0.39  |  0.22  |  NA    | NA       | NA                |   NA     |
| Smoking  | MR Egger        | 182  |0.59 (0.24-1.45)| 0.25 | 0.59           |  0.37  |        |        |          |                   |          | 
| Smoking  | Simple mode     | 182  |0.59 (0.20-1.77)| 0.35 |                |        |        |        |          |                   |          |
| Smoking  | Weighted median | 182  |0.64 (0.46-0.91)| 0.01 |                |        |        |        |          |                   |          |
| Smoking  | Weighted mode   | 182  |0.58 (0.25-1.32)| 0.20 |                |        |        |        |          |                   |          |

Here, we can see that
* **For Alcohol:**
  - There was an inverse association but not significant based on IVW random : OR=0.68, p=0.17 (there is limit significant heterogeneity showed by Q_pval = 0.06)
  - pval_intercept was not significant (p=0.94), which means that MR_Egger do not show directional pleiotropy : it is not necessary to take into account the estimate from MR-Egger
  - The weighted median and weighted mode method also showed inverse associations (not significant but very conservative)
  - The MR_Presso showed significant pleiotropy (heterogeneity, p=0.047). It found one outlier, the rs9607814. After removing it, the association was similar (OR=0.77, p=0.33) and the distortion test was not significant (p=0.71)
* **For Smoking:**
  - There was a significant inverse association showed by the IVW random : OR=0.74, p=0.01. There was no heterogeneity (q_pval=0.39, so it is possible to use as main result the IVW fix but it is recommand to use the random one)
  - pval_intercept was not significant (p=0.59), there was no directional pleiotropy
  - The weighted median and weighted mode method also showed similar associations
  - The MR-Presso did not show heterogeneity (p=0.22)


### <a id="stratMR" /> Stratified analyses and interaction

In our studies, we also implemented a new method for gxe interaction using MR. The interaction term was calculated based on the MR results from two strata (for exemple younger *versus* older). This function, to automatically write the table of results, is called "MR_inter_fun".


> MR_inter_fun(name_exposure_file, name_outcome1_file, strata1, name_outcome2_file, strata2) {}

- **name_exposure_file** is the name of the file with the associations between SNPs and exposures after clump (here: "exposures_0.001")
- **name_outcome1_file** is the name of the file with the associations between SNPs and the first outcome strata (here: "outcome_younger")
- **strata1** is the name of the strata of the outcome1 (here: "younger")
- **name_outcome2_file** is the name of the file with the associations between SNPs and the second outcome strata (here: "outcome_older")
- **strata2** is the name of the strata of the outcome2 (here: "older")
(these files are in the 0_data folder)

```
list_presso = list("Smoking", "Alcohol", "Coffee")

MR_inter_fun("exposures_0.001", "outcome_younger", "younger", "outcome_older", "older")
```

As for the MR_fun function, before running the MR_fun, you have to define the list of your exposures in the list_presso object to use them in the MR_presso function. Fill in only the exhibitions for which there are several VI (about 5 minimum). 

This function will write two tables of results: one for the first strata and one other for the second strata. They will be written in the same results file: 1_results/Results_MR_exposures_0.001.xlsx, in the "younger" and "older" page. 
The columns are the same that presented for the MR_fun, and four columns were added : beta_int, se_int, OR_int and p_int, for the estimation of the interaction. 


### <a id="robust" /> Other robust methods

Many others methods are implemented. Here is the code to estimate MR association using IVW robust, MR Lasso and MR conmix.  
This time, you have to use the MendelianRandomization Rpackage. 
The name of the function is "MR_Lasso_fun".

> MR_Lasso_fun(name_exposure_file, name_outcome_file, strata) {}

- **name_exposure_file** is the name of the file with the associations between SNPs and exposures after clump (here: "exposures_0.001")
- **name_outcome_file** is the name of the file with the associations between SNPs and outcome (here: "outcome")
- **strata** is the name of the strata of the analysis (here: "All indiv")
(these files are in the 0_data folder)

```
list_expo=list("Smoking", "Alcohol", "Coffee")
list_outcome="PD"
MR_Lasso_fun("exposures_0.001", "outcome", "All_indiv")
```

Before running the MR_fun, you have to define the list of your exposures in the list_expo object, and the list of your outcomes in the list_outcome object to use them in the MR_Lasso function. 

Running this function automatically create the Results_MRLasso_exposures_0.001.xlsx file in the 1_results folder.
This file contains several columns:

|exposure|outcome|beta|se|OR|pval|nSnps_valid|lambda|
|--------|-------|----|--|--|----|-----------|------|


### <a id="reverse" /> Reverse MR

To assess reverse causation, you also should perform a reverse MR. The methods used are the same, but the SNPs used as IVs are associated with the disease. (Actually, the disease becomes the exposure, and the exposures become the diseases)

### <a id="plots" /> MR Plots

In our study, we choose to show the scatter and funnel plot for our exposures. The MR_scatter_funnel_fun allows to do it.

> MR_scatter_funnel_fun(name_exposure_file, name_outcome_file, strata) {}

- **name_exposure_file** is the name of the file with the associations between SNPs and exposures after clump (here: "exposures_0.001")
- **name_outcome_file** is the name of the file with the associations between SNPs and outcome (here: "outcome")
- **strata** is the name of the strata of the analysis (here: "All indiv")
(these files are in the 0_data folder)

```
MR_scatter_funnel_fun("exposures_0.001","outcome","All indiv")
```

This function provides scatter and funnel plots for each exposure (.png format). These plots are saved in the 1_results folder.

Several other studs are available in the TwoSampleMR package (leave one out for exemple). Please see the documentation of the package for more details

--The END--
