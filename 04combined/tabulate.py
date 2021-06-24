# -*- coding: utf-8 -*-
# combine MAF from PLINK .frq vs R2 from INFO
import sys

r_squared = {}
mafs = {}
chromosome = {}
A1 = {}
A2 = {}

with open(sys.argv[1]) as maf_file:
    next(maf_file)
    for line in maf_file:
        name = line.split()[1]
        maf = line.split()[4]
        chro = line.split()[0]
        A1allele = line.split()[2]
        A2allele = line.split()[3]
#        print ("name:" + name + "maf:" + maf + "chromosome:" + chro + "A1:" + A1allele + "A2:"+ A2allele)
        mafs[name] = maf
        chromosome[name] = chro
        A1[name] = A1allele
        A2[name] = A2allele    


with open(sys.argv[2]) as info_file: 
    for line in info_file:
        name = line.split()[1]
        Rsq = line.split()[6]
        r_squared[name] = Rsq
        
#        print ("name: " + name + "Rsq: " + Rsq)
            

            
for i in mafs:
#    print ("i in mafs is" + i)
    print (chromosome[i] + "\t" + i + "\t" +  A1[i] + "\t" +  A2[i] + "\t" +  mafs[i] + "\t" +  r_squared[i])