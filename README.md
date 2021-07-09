# Imputation Pipeline using Minimac4 and PBWT

## Introduction

This is an in-house imputation pipeline for the genotype imputation of Singaporean Chinese, Malay and Indian datasets on the 1000G, HRC and SG10K reference panels. The pipeline uses Eagle2 for pre-phasing; Minimac4 and PBWT for genotype imputation; and  PLINK, VCFtools and BCFtools for general handling and quality control of the data. The pipeline outputs a table containing quality controlled SNPs, minor allele frequencies and imputation quality scores.

## Technologies
* Python version: 3.7.6
* Eagle2 version: 2.4.1
* Minimac3 version: 2.0.1
* Minimac4 version: 1.0.2
* PBWT version: 3.0-8c25e5c
* PLINK version: 1.90b5.3
* VCFtools version: 0.1.13
* BCFtools version: 1.9 (using htslib 1.9)
* tabix version: 0.2.5

## How to run
1. 01pre-phasing 
* Edit Parameters section in phase.sh
* PLINK should be on $PATH 
* Run with "./phase.sh".
2. 02imputation 
* Edit Parameters section in impute.sh.
* Reference panels can be converted to m3vcf and PBWT formats using Minimac3 and PBWT respectively
* Run with "./impute.sh"
3. 03postimpQC
* Edit Parameters section in postimp_qc.sh
* Run with "./postimp_qc.sh"
4. 04combined
* Edit Parameters section in combine.sh.
* VCFtools, BCFtools, bgzip and tabix should be on $PATH
* Run with "./combine.sh"
5. Final files in 04combined sorted into race and panel folders:
* "$software"_"$race"_"$panel"_chrALL_preQC.vcf.gz - .vcf.gz file without post imputation QC
* "$software"_"$race"_"$panel"_chrALL_postQC.vcf.gz - .vcf.gz file after post imputation QC 
* combined_info.txt - INFO file for all chromosomes
* maf_R2_preQC.txt - SNP list with chromosome, SNP ID, PLINK allele 1, PLINK allele 2, minor allele frequency and imputation quality score
* maf_R2_postQC.txt - SNP list with chromosome, SNP ID, PLINK allele 1, PLINK allele 2, minor allele frequency and imputation quality score after post imputation QC
* maf_R2_pre/postQC_maf1e-3.txt - same as maf_R2_pre/postQC.txt, but filtered by MAF < 0.005
