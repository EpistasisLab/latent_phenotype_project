#!/bin/bash
#BSUB -J step8.4_get_imputed_data
#BSUB -o step8.4_get_imputed_data.out
#BSUB -e step8.4_get_imputed_data.err

bsub < eid_filters/get_UKB_imputed_samples_chr1.sh
bsub < eid_filters/get_UKB_imputed_samples_chr2.sh
bsub < eid_filters/get_UKB_imputed_samples_chr3.sh
bsub < eid_filters/get_UKB_imputed_samples_chr4.sh
bsub < eid_filters/get_UKB_imputed_samples_chr5.sh
bsub < eid_filters/get_UKB_imputed_samples_chr6.sh
bsub < eid_filters/get_UKB_imputed_samples_chr7.sh
bsub < eid_filters/get_UKB_imputed_samples_chr8.sh
bsub < eid_filters/get_UKB_imputed_samples_chr9.sh
bsub < eid_filters/get_UKB_imputed_samples_chr10.sh
bsub < eid_filters/get_UKB_imputed_samples_chr11.sh
bsub < eid_filters/get_UKB_imputed_samples_chr12.sh
bsub < eid_filters/get_UKB_imputed_samples_chr13.sh
bsub < eid_filters/get_UKB_imputed_samples_chr14.sh
bsub < eid_filters/get_UKB_imputed_samples_chr15.sh
bsub < eid_filters/get_UKB_imputed_samples_chr16.sh
bsub < eid_filters/get_UKB_imputed_samples_chr17.sh
bsub < eid_filters/get_UKB_imputed_samples_chr18.sh
bsub < eid_filters/get_UKB_imputed_samples_chr19.sh
bsub < eid_filters/get_UKB_imputed_samples_chr20.sh
bsub < eid_filters/get_UKB_imputed_samples_chr21.sh
bsub < eid_filters/get_UKB_imputed_samples_chr22.sh
bsub < eid_filters/get_UKB_imputed_samples_chrMT.sh
bsub < eid_filters/get_UKB_imputed_samples_chrX.sh
bsub < eid_filters/get_UKB_imputed_samples_chrXY.sh
bsub < eid_filters/get_UKB_imputed_samples_chrY.sh
