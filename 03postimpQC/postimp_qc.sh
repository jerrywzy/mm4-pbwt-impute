#!/bin/bash
### post imputation QC

RESULTS=/hpc/home/lsiwzyj/Projects/InHouseImputation/03Imputation
POSTIMPQC=/hpc/home/lsiwzyj/Projects/InHouseImputation/04postimpQC
MERGED_REF=/hpc/home/lsiwzyj/MergedRefPanels

filterMAF="0.05"  # MAF filtering by 0.001 is done at plotting stage
filterHWE="1e-6"
filtergeno="0.05"

races=( Chinese Malay Indian )  
softwares=( minimac4 PBWT ) 
panels=( 1000G HRC SG10K SG10K_1000G SG10K_HRC ) 
CHRS=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 )  

for software in "${softwares[@]}"
do
  cd ${POSTIMPQC}
  mkdir -p "$software"
  cd "$software"
  
  for race in "${races[@]}"
  do
    cd ${POSTIMPQC}/"$software"
    mkdir -p "$race"
    cd "$race"
    
    for panel in "${panels[@]}"
    do
      cd ${POSTIMPQC}/"$software"/"$race"
      mkdir -p "$panel"
      cd "$panel"
      
      for CHR in "${CHRS[@]}"
      do
        # fix multiallelic sites and reannotate to CHR:POS:REF:ALT
        bcftools view -m2 -M2 -v snps $RESULTS/"$software"/"$race"/"$panel"/"$race"_"$panel"_chr"$CHR"_"$software".dose.vcf.gz -Oz -o temp1chr"$CHR"            
        bcftools annotate --set-id '%CHROM\:%POS\:%REF\:%ALT' temp1chr"$CHR" -Oz -o temp2chr"$CHR"
        
        # normalize indels and update imputation quality scores to IMPUTE INFO scores on PBWT files
        if [[ "$software" = "PBWT" ]]; then
          bcftools norm -d all temp2chr"$CHR" -o temp3chr"$CHR"
          bcftools +impute-info temp3chr"$CHR" -o chr"$CHR"_rmdup.dose.vcf
        else
          bcftools norm -d all temp2chr"$CHR" -o chr"$CHR"_rmdup.dose.vcf
        fi
        
        
        # convert vcf to PLINK for QC filtering
        plink --vcf chr"$CHR"_rmdup.dose.vcf --const-fid --keep-allele-order --biallelic-only strict --make-bed --out chr"$CHR"_rmdup
        
        # make info table
        if [[ "$software" = "PBWT" ]]; then
          echo "Generating INFO file from PBWT imputation"
          bcftools query -f '%CHROM %ID  %POS  %REF  %ALT %INFO/RefPanelAF %INFO/INFO %INFO/INFO\n' chr"$CHR"_rmdup.dose.vcf -o chr"$CHR".info
        else 
          echo "Generating INFO file from Minimac4 imputation"
          bcftools query -f '%CHROM %ID  %POS  %REF  %ALT %INFO/MAF %INFO/R2 %INFO/IMPUTED %INFO/TYPED %INFO/TYPED_ONLY\n' chr"$CHR"_rmdup.dose.vcf -o chr"$CHR".info
        fi
        
        # begin filtering by MAF, HWE, and call rate
        plink --bfile chr"$CHR"_rmdup --keep-allele-order --qual-scores chr"$CHR".info 7 2 --qual-threshold 0.9 --maf ${filterMAF} --hwe ${filterHWE} --geno ${filtergeno} --make-bed --out "$software"_"$race"_"$panel"_chr"$CHR"
        
        # convert to vcf again
        mkdir -p "$software"_"$race"_"$panel"_vcf_preQC   # -p means no error if existing, make if not existing
        mkdir -p "$software"_"$race"_"$panel"_vcf_postQC
        
        plink --bfile chr"$CHR"_rmdup --keep-allele-order --recode vcf --out "$software"_"$race"_"$panel"_vcf_preQC/"$software"_"$race"_"$panel"_chr"$CHR"
        plink --bfile "$software"_"$race"_"$panel"_chr"$CHR" --keep-allele-order --recode vcf --out "$software"_"$race"_"$panel"_vcf_postQC/"$software"_"$race"_"$panel"_chr"$CHR"
        
        rm temp*
      done
    done
  done
done
