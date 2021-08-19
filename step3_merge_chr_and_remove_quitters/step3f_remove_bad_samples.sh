#!/bin/bash
#BSUB -J step3f_remove_bad_samples
#BSUB -o step3f_remove_bad_samples.out
#BSUB -e step3f_remove_bad_samples.error
#BSUB -R "rusage[mem=20000MB]"
#BSUB -M 20000MB
module load plink/1.90Beta

plink --memory 15000 --bfile UKB_samples_half_filtered --remove low_quality_FIDs.txt --make-bed --out ../step4_remove_relatives/UKB_samples_filtered_with_remainder
rm UKB_samples_half_filtered.bed
rm UKB_samples_half_filtered.bim
rm UKB_samples_half_filtered.fam
rm UKB_samples_half_filtered.het
rm UKB_samples_half_filtered.imiss
rm UKB_samples_half_filtered.lmiss
rm UKB_samples_half_filtered.sexcheck
rm UKB_samples_half_filtered.log