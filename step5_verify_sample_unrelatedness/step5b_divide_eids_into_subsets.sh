#!/bin/bash
#BSUB -J step5b_divide_eids_into_subsets
#BSUB -o step5b_divide_eids_into_subsets.out
#BSUB -e step5b_divide_eids_into_subsets.error
#BSUB -R "rusage[mem=20000MB]"
#BSUB -M 20000MB

module load plink/1.90Beta6.18

plink --memory 15000 --bfile ../step4_remove_relatives/UKB_samples_unrelated --keep data_subset_content/subset1.tab --make-bed --out data_subset_content/UKB_samples_filtered_subset1
plink --memory 15000 --bfile ../step4_remove_relatives/UKB_samples_unrelated --keep data_subset_content/subset2.tab --make-bed --out data_subset_content/UKB_samples_filtered_subset2
plink --memory 15000 --bfile ../step4_remove_relatives/UKB_samples_unrelated --keep data_subset_content/subset3.tab --make-bed --out data_subset_content/UKB_samples_filtered_subset3
plink --memory 15000 --bfile ../step4_remove_relatives/UKB_samples_unrelated --keep data_subset_content/subset4.tab --make-bed --out data_subset_content/UKB_samples_filtered_subset4
plink --memory 15000 --bfile ../step4_remove_relatives/UKB_samples_unrelated --keep data_subset_content/subset5.tab --make-bed --out data_subset_content/UKB_samples_filtered_subset5
plink --memory 15000 --bfile ../step4_remove_relatives/UKB_samples_unrelated --keep data_subset_content/subset6.tab --make-bed --out data_subset_content/UKB_samples_filtered_subset6
plink --memory 15000 --bfile ../step4_remove_relatives/UKB_samples_unrelated --keep data_subset_content/subset7.tab --make-bed --out data_subset_content/UKB_samples_filtered_subset7
plink --memory 15000 --bfile ../step4_remove_relatives/UKB_samples_unrelated --keep data_subset_content/subset8.tab --make-bed --out data_subset_content/UKB_samples_filtered_subset8
plink --memory 15000 --bfile ../step4_remove_relatives/UKB_samples_unrelated --keep data_subset_content/subset9.tab --make-bed --out data_subset_content/UKB_samples_filtered_subset9
plink --memory 15000 --bfile ../step4_remove_relatives/UKB_samples_unrelated --keep data_subset_content/subset10.tab --make-bed --out data_subset_content/UKB_samples_filtered_subset10

