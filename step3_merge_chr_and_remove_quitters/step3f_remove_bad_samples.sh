#!/bin/bash
#BSUB -J step3f_remove_bad_samples
#BSUB -o step3f_remove_bad_samples.out
#BSUB -e step3f_remove_bad_samples.error
#BSUB -R "rusage[mem=20000MB]"
#BSUB -M 20000MB
module load plink/1.90Beta6.18

plink --memory 15000 --bfile UKB_samples_half_filtered --remove low_quality_FIDs.txt --make-bed --out ../step4_remove_relatives/UKB_samples_filtered_with_remainder