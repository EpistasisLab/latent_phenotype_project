#!/bin/bash
#BSUB -J step3d_remove_bad_SNPs
#BSUB -o step3d_remove_bad_SNPs.out
#BSUB -e step3d_remove_bad_SNPs.error
#BSUB -R "rusage[mem=20000MB]"
#BSUB -M 20000MB
module load plink/1.90Beta

plink --memory 15000 --bfile UKB_samples --maf 0.005 --hwe 0.000001 --geno 0.02 --make-bed --out UKB_samples_half_filtered
rm UKB_samples.bed
rm UKB_samples.bim
rm UKB_samples.fam
rm UKB_samples.frq
rm UKB_samples.hwe
rm UKB_samples.imiss
rm UKB_samples.lmiss
rm UKB_samples.log

# produces data to check for low quality samples
plink --memory 15000 --bfile UKB_samples_half_filtered --missing --out UKB_samples_half_filtered
plink --memory 15000 --bfile UKB_samples_half_filtered --het --out UKB_samples_half_filtered
plink --memory 15000 --bfile UKB_samples_half_filtered --check-sex --out UKB_samples_half_filtered
