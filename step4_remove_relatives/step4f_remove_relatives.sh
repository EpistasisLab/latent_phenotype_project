#!/bin/bash
#BSUB -J step4f_remove_relatives
#BSUB -o step4f_remove_relatives.out
#BSUB -e step4f_remove_relatives.error
#BSUB -R "rusage[mem=20000MB]"
#BSUB -M 20000MB

module load plink/1.90Beta

plink --memory 15000 --bfile UKB_samples_filtered --keep unrelated_eids.tab --make-bed --out UKB_samples_unrelated
plink --memory 15000 --bfile UKB_samples_unrelated --indep-pairwise 1000 100 0.09 --out UKB_samples_unrelated
plink --memory 15000 --bfile UKB_samples_unrelated --extract UKB_samples_unrelated.prune.in --make-bed --out ../step6_PCA/UKB_samples_unrelated_pruned
# rm UKB_samples_filtered.bed
# rm UKB_samples_filtered.bim
# rm UKB_samples_filtered.fam
# rm UKB_samples_filtered.log


