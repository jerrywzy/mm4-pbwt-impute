#!/bin/sh
### impute with Minimac4 and PWBT

################
## PARAMETERS ##
################
MINIMAC4=/hpc/home/lsiwzyj/programs/Minimac4/release-build/
PBWT=/hpc/home/lsiwzyj/programs/pbwt

working=/hpc/home/lsiwzyj/Projects/InHouseImputation   # directory with main folders (01pre-phasing, 02imputation etc)
TH_REF=/hpc/home/lsiwzyj/1000GData/reference_panels    # directory with 1000G reference panels in m3vcf AND vcf format
TH_REF_PBWT=${TH_REF}/PBWT_convert    # 1000G panels in PBWT format
SG10K_REF=/hpc/home/lsiwzyj/SG10KData  # SG10K panels in VCF format
SG10K_REF_PBWT=${SG10K_REF}/pbwt   # SG10K panels in PBWT format
SG10K_REF_MM=${SG10K_REF}/m3vcf   # SG10K panels in m3vcf format
HRC_REF=/hpc/home/lsiwzyj/HRC1.1Data   # HRC panels root folder
HRC_REF_PBWT=${HRC_REF}/pbwt   # HRC panels in PBWT format
HRC_REF_MM=${HRC_REF}/m3vcf    # HRC panels in m3vcf format
PHASING=${working}/02pre-phasing
IMP=${working}/03Imputation
MERGED_REF=/hpc/home/lsiwzyj/MergedRefPanels    # merged reference panels directory for 1000G+SG10K and HRC+SG10K

softwares=( minimac4 PBWT )   
races=( Chinese Malay Indian )  
panels=( 1000G HRC SG10K SG10K_1000G SG10K_HRC )
CHRS=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 )  


################
## IMPUTATION ##
################
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
            $PBWT/pbwt -checkpoint 100000 -readVcfGT $PHASING/"$race"/"$race"_"$panel"_chr"$CHR"_phased.vcf -referenceImpute $TH_REF_PBWT/ALL.chr"$CHR".phase3_v5.shapeit2_mvncall_integrated.noSingleton.genotypes -writeVcfGz "$race"_"$panel"_chr"$CHR"_PBWT.dose.vcf.gz  
          
          elif [[ "$panel" = "SG10K" ]]; then
            echo "Panel: SG10K"
            $PBWT/pbwt -checkpoint 100000 -readVcfGT $PHASING/"$race"/"$race"_HRC_chr"$CHR"_phased.vcf -referenceImpute $SG10K_REF_PBWT/chr"$CHR" -writeVcfGz "$race"_"$panel"_chr"$CHR"_PBWT.dose.vcf.gz
          
          elif [[ "$panel" = "HRC" ]]; then
            echo "Panel: HRC"
            $PBWT/pbwt -checkpoint 100000 -readVcfGT $PHASING/"$race"/"$race"_HRC_chr"$CHR"_phased.vcf -referenceImpute $HRC_REF_PBWT/HRC.r1-1.EGA.GRCh37.chr"$CHR".haplotypes -writeVcfGz "$race"_"$panel"_chr"$CHR"_PBWT.dose.vcf.gz      
          
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
            $MINIMAC4/minimac4 --refHaps $SG10K_REF_MM/chr"$CHR".m3vcf.gz \
            --haps $PHASING/"$race"/"$race"_HRC_chr"$CHR"_phased.vcf \
            --prefix "$race"_"$panel"_chr"$CHR"_minimac4
          
          elif [[ "$panel" = "HRC" ]]; then
            echo "Panel: HRC"
            $MINIMAC4/minimac4 --refHaps $HRC_REF_MM/HRC.r1-1.EGA.GRCh37.chr"$CHR".haplotypes.m3vcf.gz \
            --haps $PHASING/"$race"/"$race"_HRC_chr"$CHR"_phased.vcf \
            --prefix "$race"_"$panel"_chr"$CHR"_minimac4
                    
          elif [[ "$panel" = "SG10K_1000G" ]]; then             
            echo "Panel: SG10K_1000G"
            $MINIMAC4/minimac4 --refHaps $MERGED_REF/SG10K_1000G/merge_3/m3vcf/SG10K_1000G_chr"$CHR"_IDfixed_overlap_merged_geno_recode.m3vcf.gz \
            --haps $PHASING/"$race"/"$race"_HRC_chr"$CHR"_phased.vcf \
            --prefix "$race"_"$panel"_chr"$CHR"_minimac4
          
          elif [[ "$panel" = "SG10K_HRC" ]]; then             
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

