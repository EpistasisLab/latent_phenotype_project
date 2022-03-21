#!/bin/bash
#BSUB -J step10b_get_significant_SNPs
#BSUB -o step10b_get_significant_SNPs.out
#BSUB -e step10b_get_significant_SNPs.err
#BSUB -R "rusage[mem=20000MB]"
#BSUB -M 20000MB

module load plink/1.90Beta

IN_PREFIX="../step8_get_imputed_ukb_samples/filtered_output/UKB_samples_chr"
for NAME in "_alcohol" "_smoking" "_gender" "_exercise" "_age" "_BMI" "_pol"
do

    DIR="significant_SNPs_plink_files"$NAME
    if test ! -d $DIR; then
        mkdir "$DIR"
    fi

    OUT_PREFIX=$DIR"/significant_SNPs_chr"

    for i in {1..22} "MT" "X" "XY" "Y"
    do
       plink --memory 15000 --bfile $IN_PREFIX$i --extract significant_rsIDs$NAME.txt --make-bed --out $OUT_PREFIX$i
    done

done