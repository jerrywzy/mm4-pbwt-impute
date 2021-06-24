# -*- coding: utf-8 -*-
### combine minor allele frequencies from PLINK .frq files with corresponding R2 scores from INFO files
import sys

r_squared = {}
mafs = {}
chromosome = {}
A1 = {}
A2 = {}

# read in .frq files
with open(sys.argv[1]) as maf_file:
    next(maf_file)
    for line in maf_file:
        name = line.split()[1]
        maf = line.split()[4]
        chro = line.split()[0]
        A1allele = line.split()[2]
        A2allele = line.split()[3]
        mafs[name] = maf
        chromosome[name] = chro
        A1[name] = A1allele
        A2[name] = A2allele    

# read in INFO file
with open(sys.argv[2]) as info_file: 
    for line in info_file:
        name = line.split()[1]
        Rsq = line.split()[6]
        r_squared[name] = Rsq
        
# print final table     
for i in mafs:
    print (chromosome[i] + "\t" + i + "\t" +  A1[i] + "\t" +  A2[i] + "\t" +  mafs[i] + "\t" +  r_squared[i])