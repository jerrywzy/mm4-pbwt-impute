#!/bin/sh

# impute with mm4 and PWBT


MINIMAC4=/hpc/home/lsiwzyj/programs/Minimac4/release-build/
PBWT=/hpc/home/lsiwzyj/programs/pbwt

TH_REF=/hpc/home/lsiwzyj/1000GData/reference_panels
SG10K_REF=/hpc/home/lsiwzyj/SG10KData
HRC_REF=/hpc/home/lsiwzyj/HRC1.1Data
PHASING=/hpc/home/lsiwzyj/Projects/InHouseImputation/02pre-phasing
IMP=/hpc/home/lsiwzyj/Projects/InHouseImputation/03Imputation
MERGED_REF=/hpc/home/lsiwzyj/MergedRefPanels

softwares=( minimac4 PBWT ) #  
races=( Chinese Malay Indian )  
panels=( SG10K_HRC ) # SG10K_1000G  taken out 1000G (done) and SG10K (done) HRC (done) 
CHRS=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 )  ## 1..22
           
for software in "${softwares[@]}"
do

  cd ${IMP}
  mkdir -p "$software"
  cd "$software"

  for race in "${races[@]}"
  do
  
    cd ${IMP}/"$software"
    mkdir -p "$race"
    cd "$race"
    
    for panel in "${panels[@]}"
    do
    
      cd ${IMP}/"$software"/"$race"
      mkdir -p "$panel"
      cd "$panel"
      
      for CHR in "${CHRS[@]}"
      do

        if [[ "$software" = "PBWT" ]]; then
          echo "Running PBWT imputation"
          
          if [[ "$panel" = "1000G" ]]; then
            echo "Panel: 1000G"
            $PBWT/pbwt -checkpoint 100000 -readVcfGT $PHASING/"$race"/"$race"_"$panel"_chr"$CHR"_phased.vcf -referenceImpute $TH_REF/PBWT_convert/ALL.chr"$CHR".phase3_v5.shapeit2_mvncall_integrated.noSingleton.genotypes -writeVcfGz "$race"_"$panel"_chr"$CHR"_PBWT.dose.vcf.gz
            
          elif [[ "$panel" = "SG10K" ]]; then
            echo "Panel: SG10K"
            $PBWT/pbwt -checkpoint 100000 -readVcfGT $PHASING/"$race"/"$race"_HRC_chr"$CHR"_phased.vcf -referenceImpute $SG10K_REF/pbwt/chr"$CHR" -writeVcfGz "$race"_"$panel"_chr"$CHR"_PBWT.dose.vcf.gz
            # note: SG10K imputation here uses HRC phased file from 6.5 further exclusion, as there is no equivalent file from Mich server where differences in SNPs are excluded. I chose HRC as it is has the least SNPs excluded in 6.5
          elif [[ "$panel" = "HRC" ]]; then
            echo "Panel: HRC"
            $PBWT/pbwt -checkpoint 100000 -readVcfGT $PHASING/"$race"/"$race"_HRC_chr"$CHR"_phased.vcf -referenceImpute $HRC_REF/pbwt/HRC.r1-1.EGA.GRCh37.chr"$CHR".haplotypes -writeVcfGz "$race"_"$panel"_chr"$CHR"_PBWT.dose.vcf.gz
            
          elif [[ "$panel" = "1000G_noEUR" ]]; then
            echo "Panel: 1000G_noEUR"
            $PBWT/pbwt -checkpoint 100000 -readVcfGT $PHASING/"$race"/"$race"_HRC_chr"$CHR"_phased.vcf -referenceImpute /hpc/home/lsiwzyj/1000GData/subsets/excluded_EUR/pbwt/ALL.chr"$CHR".phase3_v5.shapeit2_mvncall_integrated.noSingleton.genotypes.noEUR -writeVcfGz "$race"_"$panel"_chr"$CHR"_PBWT.dose.vcf.gz
            
          elif [[ "$panel" = "SG10K_1000G" ]]; then
            echo "Panel: SG10K_1000G"
            $PBWT/pbwt -checkpoint 100000 -readVcfGT $PHASING/"$race"/"$race"_HRC_chr"$CHR"_phased.vcf -referenceImpute $MERGED_REF/SG10K_1000G/merge_3/pbwt/SG10K_1000G_chr"$CHR"_IDfixed_overlap_merged_geno_recode -writeVcfGz "$race"_"$panel"_chr"$CHR"_PBWT.dose.vcf.gz
            
          elif [[ "$panel" = "SG10K_HRC" ]]; then
            echo "Panel: SG10K_HRC"
            $PBWT/pbwt -checkpoint 100000 -readVcfGT $PHASING/"$race"/"$race"_HRC_chr"$CHR"_phased.vcf -referenceImpute $MERGED_REF/SG10K_HRC/merge_3/pbwt/SG10K_HRC_chr"$CHR"_IDfixed_overlap_merged_geno_recode -writeVcfGz "$race"_"$panel"_chr"$CHR"_PBWT.dose.vcf.gz
            
          fi
        elif [[ "$software" = "minimac4" ]]; then 
          echo "Running Minimac4 imputation"
          
          if [[ "$panel" = "1000G" ]]; then
            $MINIMAC4/minimac4 --refHaps $TH_REF/"$CHR".1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz \
            --haps $PHASING/"$race"/"$race"_"$panel"_chr"$CHR"_phased.vcf \
            --prefix "$race"_"$panel"_chr"$CHR"_minimac4
          elif [[ "$panel" = "SG10K" ]]; then
            echo "Panel: SG10K"
            $MINIMAC4/minimac4 --refHaps $SG10K_REF/m3vcf/chr"$CHR".m3vcf.gz \
            --haps $PHASING/"$race"/"$race"_HRC_chr"$CHR"_phased.vcf \
            --prefix "$race"_"$panel"_chr"$CHR"_minimac4
            
          elif [[ "$panel" = "HRC" ]]; then
            echo "Panel: HRC"
            $MINIMAC4/minimac4 --refHaps $HRC_REF/m3vcf/HRC.r1-1.EGA.GRCh37.chr"$CHR".haplotypes.m3vcf.gz \
            --haps $PHASING/"$race"/"$race"_HRC_chr"$CHR"_phased.vcf \
            --prefix "$race"_"$panel"_chr"$CHR"_minimac4
            
          elif [[ "$panel" = "1000G_noEUR" ]]; then             # NOT DONE
            echo "Panel: 1000G_noEUR"
            $MINIMAC4/minimac4 --refHaps $xxx.m3vcf.gz \
            --haps $PHASING/"$race"/"$race"_HRC_chr"$CHR"_phased.vcf \
            --prefix "$race"_"$panel"_chr"$CHR"_minimac4
            
          elif [[ "$panel" = "SG10K_1000G" ]]; then             #  
            echo "Panel: SG10K_1000G"
            $MINIMAC4/minimac4 --refHaps $MERGED_REF/SG10K_1000G/merge_3/m3vcf/SG10K_1000G_chr"$CHR"_IDfixed_overlap_merged_geno_recode.m3vcf.gz \
            --haps $PHASING/"$race"/"$race"_HRC_chr"$CHR"_phased.vcf \
            --prefix "$race"_"$panel"_chr"$CHR"_minimac4
            
          elif [[ "$panel" = "SG10K_HRC" ]]; then             #  
            echo "Panel: SG10K_HRC"
            $MINIMAC4/minimac4 --refHaps $MERGED_REF/SG10K_HRC/merge_3/m3vcf/SG10K_HRC_chr"$CHR"_IDfixed_overlap_merged_geno_recode.m3vcf.gz \
            --haps $PHASING/"$race"/"$race"_HRC_chr"$CHR"_phased.vcf \
            --prefix "$race"_"$panel"_chr"$CHR"_minimac4
          
          fi
        fi
      
      done
  
    done

  done
     
done

