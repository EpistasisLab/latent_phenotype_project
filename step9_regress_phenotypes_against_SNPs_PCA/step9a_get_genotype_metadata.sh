#!/bin/bash
#BSUB -J step9a_make_genotype_metadata
#BSUB -o step9a_make_genotype_metadata.out
#BSUB -e step9a_make_genotype_metadata.err
#BSUB -R "rusage[mem=20000MB]"
#BSUB -M 20000MB

module load plink/1.90Beta6.18

DIR="genotype_metadata"
if test ! -d $DIR; then
    mkdir "$DIR"
fi

IN_PREFIX="../step8_get_imputed_ukb_samples/filtered_output/UKB_samples_chr"
OUT_PREFIX="genotype_metadata/genotype_metadata_chr"
for i in {1..22} "MT" "X" "XY" "Y"
do
   plink --memory 15000 --bfile $IN_PREFIX$i --freqx  --out $OUT_PREFIX$i
   plink --memory 15000 --bfile $IN_PREFIX$i --keep-allele-order --freq  --out $OUT_PREFIX$i
done

# plink --memory 15000 --bfile ../step8_get_imputed_ukb_samples/filtered_output/UKB_samples_chrY --keep-allele-order--freqx --out genotype_metadata/chrY