#!/bin/sh
### phase input files with Eagle2

################
## PARAMETERS ##
################
working=/hpc/home/lsiwzyj/Projects/InHouseImputation/   # directory with main folders (01pre-phasing, 02imputation etc)
EAGLE=/hpc/home/lsiwzyj/programs/Eagle_v2.4.1  # directory containing the Eagle executable
TARGET=/hpc/home/lsiwzyj/iOmics_QC/Genomic-Feb2020/Genotyping/Omni2.5    # directory containing QC-ed VCF input files
TABLES=$EAGLE/tables    # directory with genetic_map_hg19_withX.txt.gz, included in the Eagle tar file in /tables folder
PHASING=${working}/02pre-phasing

###########
## Phase ##
###########
# Chinese
cd ${PHASING}
mkdir -p Chinese
cd Chinese

for CHR in {1..22}
do
  $EAGLE/eagle --vcf=$TARGET/Chinese_HRC/06.5further_exclusion/1000G/iOmics_Chinese_110samples_autosomes-updated-chr"$CHR"-exclusion.recode.vcf --chrom="$CHR" --geneticMapFile=$TABLES/genetic_map_hg19_withX.txt.gz --vcfOutFormat=v --outPrefix=Chinese_1000G_chr"$CHR"_phased
  
  $EAGLE/eagle --vcf=$TARGET/Chinese_HRC/06.5further_exclusion/HRC/iOmics_Chinese_110samples_autosomes-updated-chr"$CHR"-exclusion.recode.vcf --chrom="$CHR" --geneticMapFile=$TABLES/genetic_map_hg19_withX.txt.gz --vcfOutFormat=v --outPrefix=Chinese_HRC_chr"$CHR"_phased
done

# Malay
cd ${PHASING}
mkdir -p Malay
cd Malay

for CHR in {1..22}
do
  $EAGLE/eagle --vcf=$TARGET/Malay_HRC/06.5further_exclusion/1000G/iOmics_Malay_108samples_autosomes-updated-chr"$CHR"-exclusion.recode.vcf --chrom="$CHR" --geneticMapFile=$TABLES/genetic_map_hg19_withX.txt.gz --vcfOutFormat=v --outPrefix=Malay_1000G_chr"$CHR"_phased
  
  $EAGLE/eagle --vcf=$TARGET/Malay_HRC/06.5further_exclusion/HRC/iOmics_Malay_108samples_autosomes-updated-chr"$CHR"-exclusion.recode.vcf --chrom="$CHR" --geneticMapFile=$TABLES/genetic_map_hg19_withX.txt.gz --vcfOutFormat=v --outPrefix=Malay_HRC_chr"$CHR"_phased
done

# Indian
cd ${PHASING}
mkdir -p Indian
cd Indian

for CHR in {1..22}
do
  $EAGLE/eagle --vcf=$TARGET/Indian_HRC/06.5further_exclusion/1000G/iOmics_Indian_105samples_autosomes-updated-chr"$CHR"-exclusion.recode.vcf --chrom="$CHR" --geneticMapFile=$TABLES/genetic_map_hg19_withX.txt.gz --vcfOutFormat=v --outPrefix=Indian_1000G_chr"$CHR"_phased
  
  $EAGLE/eagle --vcf=$TARGET/Indian_HRC/06.5further_exclusion/HRC/iOmics_Indian_105samples_autosomes-updated-chr"$CHR"-exclusion.recode.vcf --chrom="$CHR" --geneticMapFile=$TABLES/genetic_map_hg19_withX.txt.gz --vcfOutFormat=v --outPrefix=Indian_HRC_chr"$CHR"_phased
done