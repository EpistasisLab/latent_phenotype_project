#!/bin/bash
#BSUB -J step6b_run_PCA
#BSUB -o step6b_run_PCA.out
#BSUB -e step6b_run_PCA.error
#BSUB -R "rusage[mem=150000MB]"
#BSUB -M 150000MB
module load eigensoft/6.0.1

/appl/eigensoft-6.0.1/bin/smartpca -p step0_run_PCA.par
