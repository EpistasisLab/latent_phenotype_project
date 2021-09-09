#!/bin/bash
#BSUB -J step10b_get_significant_SNPs
#BSUB -o step10b_get_significant_SNPs.out
#BSUB -e step10b_get_significant_SNPs.err
#BSUB -R "rusage[mem=20000MB]"
#BSUB -M 20000MB

module load plink/1.90Beta

DIR="significant_SNPs_plink_files"
if test ! -d $DIR; then
    mkdir "$DIR"
fi

IN_PREFIX="../step8_get_imputed_ukb_samples/filtered_output/UKB_samples_chr"
OUT_PREFIX="significant_SNPs_plink_files/significant_SNPs_chr"
for i in {1..22} "MT" "X" "XY" "Y"
do
   plink --memory 15000 --bfile $IN_PREFIX$i --extract significant_rsIDs.txt --make-bed --out $OUT_PREFIX$i
done