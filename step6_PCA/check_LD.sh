#!/bin/bash
#BSUB -J check_LD
#BSUB -o check_LD.out
#BSUB -e check_LD.error
#BSUB -R "rusage[mem=20000MB]"
#BSUB -M 20000MB
module load plink/1.90Beta

plink --memory 15000 --bfile UKB_samples_unrelated_pruned --r --keep-allele-order --ld-window-kb 200 --ld-window 1000000 --with-freqs --ld-window-r2 0 --out UKB_samples_unrelated_pruned
