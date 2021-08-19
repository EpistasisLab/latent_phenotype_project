#!/bin/bash
#BSUB -J step3b_merge_datasets
#BSUB -o step3b_merge_datasets.out
#BSUB -e step3b_merge_datasets.error
#BSUB -R "rusage[mem=20000MB]"
#BSUB -M 20000MB
module load plink/1.90Beta

# makes main genotype files
plink --memory 15000 --bfile ../step2_get_UKB_samples/filtered_output/UKB_samples_chr1 --merge-list step3.0_output_file_names.txt --remove people_who_quit.txt --make-bed --out UKB_samples

# produces data to check for low quality SNPs. 
plink --memory 15000 --bfile UKB_samples --freq --out UKB_samples
plink --memory 15000 --bfile UKB_samples --missing --out UKB_samples
plink --memory 15000 --bfile UKB_samples --hardy --out UKB_samples





