#!/bin/bash

VCFTOOLS=/hpc/home/lsiwzyj/perl5/lib64/perl5/vcftools-vcftools-cb8e254/src/perl
BGZIP=/hpc/home/lsiwzyj/programs/bin

POSTIMPQC=/hpc/home/lsiwzyj/Projects/InHouseImputation/04postimpQC
COMBINED=/hpc/home/lsiwzyj/Projects/InHouseImputation/05combined

races=( Chinese Malay Indian )  # Chinese Malay Indian
softwares=( minimac4 PBWT ) #    ### 
panels=( SG10K_HRC ) # SG10K_1000G   1000G HRC SG10K

## settle directories
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
      
      # concat vcfs for pre and post imp qc
      bcftools concat -o concat_preQC.vcf $POSTIMPQC/"$software"/"$race"/"$panel"/"$software"_"$race"_"$panel"_vcf_preQC/"$software"_"$race"_"$panel"_chr*.vcf           
      bcftools concat -o concat_postQC.vcf $POSTIMPQC/"$software"/"$race"/"$panel"/"$software"_"$race"_"$panel"_vcf_postQC/"$software"_"$race"_"$panel"_chr*.vcf

      ## sort and zip for final file keeping (this vcf.gz is not used further after this point)
      vcf-sort concat_preQC.vcf > "$software"_"$race"_"$panel"_chrALL_preQC.vcf
      bgzip -c "$software"_"$race"_"$panel"_chrALL_preQC.vcf > "$software"_"$race"_"$panel"_chrALL_preQC.vcf.gz
      tabix -p vcf "$software"_"$race"_"$panel"_chrALL_preQC.vcf.gz
      
      vcf-sort concat_postQC.vcf > "$software"_"$race"_"$panel"_chrALL_postQC.vcf
      bgzip -c "$software"_"$race"_"$panel"_chrALL_postQC.vcf > "$software"_"$race"_"$panel"_chrALL_postQC.vcf.gz
      tabix -p vcf "$software"_"$race"_"$panel"_chrALL_postQC.vcf.gz

      ## get combined .info table, this is ALL sites
      cat $POSTIMPQC/"$software"/"$race"/"$panel"/*info > combined_info.txt

      # make plink .frq file from combined VCF for pre and post QC
      plink --vcf "$software"_"$race"_"$panel"_chrALL_preQC.vcf --const-fid --freq --out freqs_chrALL_preQC
      plink --vcf "$software"_"$race"_"$panel"_chrALL_postQC.vcf --const-fid --freq --out freqs_chrALL_postQC

      # python script to create maf v R2 table for plotting
      python ${COMBINED}/tabulate.py ${COMBINED}/"$software"/"$race"/"$panel"/freqs_chrALL_preQC.frq ${COMBINED}/"$software"/"$race"/"$panel"/combined_info.txt > maf_R2_preQC.txt
      python ${COMBINED}/tabulate.py ${COMBINED}/"$software"/"$race"/"$panel"/freqs_chrALL_postQC.frq ${COMBINED}/"$software"/"$race"/"$panel"/combined_info.txt > maf_R2_postQC.txt
      
      
      # smaller filtered by maf table
      awk '$5 > 0.001' maf_R2_preQC.txt > maf_R2_preQC_maf1e-3.txt
      awk '$5 > 0.001' maf_R2_postQC.txt > maf_R2_postQC_maf1e-3.txt
      
    done
    
  done
  
done
