#!/bin/bash
### combine imputed files by chromosome, and create files filtered by MAF < 0.001

################
## PARAMETERS ##
################
working=/hpc/home/lsiwzyj/Projects/InHouseImputation/   # directory with main folders (01pre-phasing, 02imputation etc)
POSTIMPQC=${working}/04postimpQC
COMBINED=${working}/05combined

races=( Chinese Malay Indian ) 
softwares=( minimac4 PBWT ) 
panels=( 1000G HRC SG10K SG10K_1000G SG10K_HRC )   

#############
## COMBINE ##
#############
for software in "${softwares[@]}"
do
  cd ${COMBINED}
  mkdir -p "$software"
  cd "$software"
  
  for race in "${races[@]}"
  do
    cd ${COMBINED}/"$software"
    mkdir -p "$race"
    cd "$race"
    
    for panel in "${panels[@]}"
    do
      cd ${COMBINED}/"$software"/"$race"
      mkdir -p "$panel"
      cd "$panel"
      
      # concatenate VCFs from pre and post imputation QC
      bcftools concat -o concat_preQC.vcf $POSTIMPQC/"$software"/"$race"/"$panel"/"$software"_"$race"_"$panel"_vcf_preQC/"$software"_"$race"_"$panel"_chr*.vcf           
      bcftools concat -o concat_postQC.vcf $POSTIMPQC/"$software"/"$race"/"$panel"/"$software"_"$race"_"$panel"_vcf_postQC/"$software"_"$race"_"$panel"_chr*.vcf

      # sort and bgzip for final file keeping
      vcf-sort concat_preQC.vcf > "$software"_"$race"_"$panel"_chrALL_preQC.vcf
      bgzip -c "$software"_"$race"_"$panel"_chrALL_preQC.vcf > "$software"_"$race"_"$panel"_chrALL_preQC.vcf.gz
      tabix -p vcf "$software"_"$race"_"$panel"_chrALL_preQC.vcf.gz
      
      vcf-sort concat_postQC.vcf > "$software"_"$race"_"$panel"_chrALL_postQC.vcf
      bgzip -c "$software"_"$race"_"$panel"_chrALL_postQC.vcf > "$software"_"$race"_"$panel"_chrALL_postQC.vcf.gz
      tabix -p vcf "$software"_"$race"_"$panel"_chrALL_postQC.vcf.gz

      # get concatenated .info table from all chromosomes
      cat $POSTIMPQC/"$software"/"$race"/"$panel"/*info > combined_info.txt

      # create plink .frq file from combined VCF for pre and post QC
      plink --vcf "$software"_"$race"_"$panel"_chrALL_preQC.vcf --const-fid --freq --out freqs_chrALL_preQC
      plink --vcf "$software"_"$race"_"$panel"_chrALL_postQC.vcf --const-fid --freq --out freqs_chrALL_postQC

      # run tabulate.py python script to create maf v R2 table for plotting
      python ${COMBINED}/tabulate.py ${COMBINED}/"$software"/"$race"/"$panel"/freqs_chrALL_preQC.frq ${COMBINED}/"$software"/"$race"/"$panel"/combined_info.txt > maf_R2_preQC.txt
      python ${COMBINED}/tabulate.py ${COMBINED}/"$software"/"$race"/"$panel"/freqs_chrALL_postQC.frq ${COMBINED}/"$software"/"$race"/"$panel"/combined_info.txt > maf_R2_postQC.txt
      
      # create smaller filtered by maf table using awk
      awk '$5 > 0.001' maf_R2_preQC.txt > maf_R2_preQC_maf1e-3.txt
      awk '$5 > 0.001' maf_R2_postQC.txt > maf_R2_postQC_maf1e-3.txt
    done
  done
done
