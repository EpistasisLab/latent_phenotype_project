#!/bin/bash

rsIDs=("rs1906609" "rs7857118" "rs12627426" "rs73839819" "rs2234962" "rs12138073" "rs73188900")

# Array of outfiles
outfiles=("figure3b.txt" "figure3b_TRACE_PCA.txt" "figure3b_TRACE_logistic_PCA.txt" "figure3b_TRACE_NN.txt")

# Array of paths
paths=(
  "/home/greggj/pleiotropy_and_GxE/step9_regress_phenotypes_against_SNPs_logistic_PCA/binary_HF_QTL_output/*"
  "/home/greggj/pleiotropy_and_GxE/step9_regress_phenotypes_against_SNPs_PCA/QTL_output_smoking/*"
  "/home/greggj/pleiotropy_and_GxE/step9_regress_phenotypes_against_SNPs_logistic_PCA/QTL_output_smoking/*"
  "/home/greggj/pleiotropy_and_GxE/step9_regress_phenotypes_against_SNPs_NN/QTL_output_smoking/*"
)

# Loop over the indices of the outfiles array
for index in ${!outfiles[*]}
do
    outfile=${outfiles[$index]}
    path=${paths[$index]}

    if [ -f $outfile ] ; then
        rm $outfile
    fi

    for rsID in ${rsIDs[@]}
    do
        # Find all files in the directory
        for file in $path
        do
            # Use grep with word boundaries (\b) to find exact rsID matches
            grep -P "\b${rsID}\b" $file >> $outfile
        done
    done
done



