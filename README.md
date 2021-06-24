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

## Setup
