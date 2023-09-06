"# latent_phenotype_project" 

IMPORTANT: File names are listed and described in the order that they are supposed to be run. 

IMPORTANT: Any file name containing "step0" is not meant to be run directly. Rather, some other file in the sequence calls that file to run when necessary.

## Directory: step1_get_phenotypes_complex

 - `step1a_get_UKB_phenotypes.py`: imports data fields from our our UKB data copies. Selects the first measured value for each field when multiple measurements are taken over time. Modified pack years of smoking to be 0 for nonsmokers and missing otherwise. Computes each individual's total annual alcohol consumption. Converts categorical variables into binary variables. Replaces certain missing values with 0 when appropriate.

 - `step1a_library.py`: contains functions used by step1a_get_UKB_phenotypes.py 

## Directory: step1_get_phenotypes_simple

 - `step1a_import_HF_ICD_codes_and_covariates.py`: Imports more fields. Binarizes some and renames others.

 - `step1b_process_HF_ICD_codes_and_covariates.py`: Binarizes the ICD code statuses for each patient (1 for present, 0 for absent). Creates the all cause heart failure variable (referred to as "any_HF"). Removes ICD code columns not significantly correlated to all cause heart failure.

## Directory: step2_get_UKB_samples

 - `step2a_make_UKB_sample_getters.py`: makes bash scripts that read in non-imputed genotype data for each individual.

 - `step2b_get_UKB_samples.sh`: runs all bash scripts generated in the previous step

## Directory: step3_merge_chr_and_remove_quitters

 - `step3a_get_people_who_quit.py`: creates a list of eids corresponding to individuals who disallowed the UK biobank from continuing to use their data.

 - `step3b_merge_datasets.sh`: merges the chromosomes' plink files from step 2 into a single file. Removes individuals who quit according to the list from step3a.

 - `step3c_find_SNP_cutoffs.py`: generates figures related to SNPs' MAFs, individual missingness tates, and hardy weinberg equilibrium p values to help us determine removal thresholds. 

 - `step3d_remove_bad_SNPs.sh`: removes all SNPs exceeding our selected MAF, hardy weinberg equilibrium p value, or SNP missingness rate thresholds.

 - `step3e_find_sample_cutoffs.py`: generates figures related to individuals' heterozygosity, average SNP missingness rates, and X chromosome heterozygosity F scores.  

 - `step3f_remove_bad_samples.sh`: removes all individuals exceeding our selected heterozygosity, X chromosome heterozygosity F score, or average SNP missingness thresholds. 

## Directory: step4_remove_relatives

 - `step4a_divide_eids_into_subsets_prep.py`: divides the eids remaining after step 3 into 10 equally sized subsets. 

 - `step4b_divide_eids_into_subsets.sh`: Divides the SNPs into ten subsets as determined by the ten eid subsets in the previous step. The remainder eids are removed.

 - `step4c_make_king_subjobs.py`: Based on the ten new SNP subset plink files, KING subjob bash scripts are generated with this python file. 

 - `step4d_run_king_subjobs.sh`: This file runs the KING subjob bash scripts that the previous step generated. 

 - `step4e_get_unrelated_eids.py`: This file takes all pairs of individuals that are third degree relatives or higher and prunes individuals from the dataset such that, in order of importance, 1) no pair of individuals in the remaining dataset are higher than fourth degree relatives, 2) the number of all cause heart failure cases is maximized, and 3) the number of total individuals is maximized. This means that fewer total individuals are accepted when doing so increases the number of cases. 

 - `step4f_remove_relatives.sh`: removes the individuals that were selected in the previous step to ensure that GWAS is conducted on unrelated individuals. Also creates a pruned version of this dataset with a subset of SNPs in low LD with one-another to be used in computing the SNPs' principal components to adjust for population stratification.  

## Directory: step5_verify_sample_unrelatedness

 - `step5a_divide_eids_into_subsets_prep.py`: For the purpose of validating step 4, divides the eids remaining after step 3 into 10 equally sized subsets. 

 - `step5b_divide_eids_into_subsets.sh`: For the purpose of validating step 4, divides the SNPs into ten subsets as determined by the ten eid subsets in the previous step.

 - `step5c_make_king_subjobs.py`: For the purpose of validating step 4, based on the ten new SNP subset plink files, KING subjob bash scripts are generated with this python file. 

 - `step5d_run_king_subjobs.sh`: For the purpose of validating step 4, this file runs the KING subjob bash scripts that the previous step generated.  

 - `step5e_confirm_unrelatedness.py`: For the purpose of validating step 4, this file confirms that the output from step5d_run_king_subjobs.sh has identified no related individuals from the output of step 4. 

## Directory: step6_PCA

 - `check_LD.sh`: Computes the LD between SNP pairs from the pruned SNP set generated by step4f_remove_relatives.sh

 - `check_LD.py`: Confirms that no SNP pairs' R^2 value output by check_LD.sh exceeds the threshold set in step4f_remove_relatives.sh

 - `step0_run_PCA.par`: specifies parameters for eigensoft to compute principal components 
 
 - `step6a_make_PCA_documents.py`: generates files for eigensoft to compute principal components

 - `step6b_run_PCA.sh`: tells eigensoft to compute principal components

## Directory: step7_adjust_HF_for_covariates_logistic_PCA

 - `step7a_get_HF_ICD_codes_unrelated.py`: Imports phenotypes and relevant UK biobank fields corresponding to the eids of unrelated individuals determined in step 4.

 - `step7b_logPCA_transform_the_data.py`: DO NOT RUN DIRECTLY. This file applies logistic PCA dimensionality reduction to 311 ICD10 codes and all cause heart failure. The number of latent features to generate is specified via argparse with parameter k. Also generates the 16th latent phenotype defined by the all-cause heart failure residuals. Adjusts all latent phenotypes by the principal components from step 6.

 - `step7b_logPCA_transform_the_data.sh`: This file runs step7b_logPCA_transform_the_data.py to create latent features for k = 1, 2, ..., 20. We have ended up using k = 15 in the paper. 

 - `step7c_impute_missing_values.py`: DO NOT RUN DIRECTLY. This file applies MICE imputation to environmental factors' missing values. Features are filtered based on their strongest correlation with the environmental factors that missing values are being imputed for. A threshold for this correlation is specified via argparse parameter "nn". Also simulates different types of missingness to measure the imputation accuracy. 

 - `step7c_impute_missing_values.sh`: This file runs step7c_impute_missing_values.py to create latent features for nn = 0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, and 0.1. We ended up using 0.05 in the paper. 

 - `step7d_get_imputed_values_and_trasformed_variables.py`: essentially runs like step7c_impute_missing_values.py. The argparse nn parameter is not included, but instead set at 0.05. Also missingness is not simulated this time. Rather, missing values are imputed, and environmental factors with imputed missing values are produced. 

 - `step7d_get_imputed_values_and_trasformed_variables.sh`: Runs step7d_get_imputed_values_and_trasformed_variables.py. You could just run step7d_get_imputed_values_and_trasformed_variables.py directly, but it's runtime is such that submitting a job is less prone to disruption. 

 - `step7e_get_PCs_effs.py`: Produces files with genetic principle components and beta coefficients for logistic regression between all cause heart failure and the genetic principle components. Used for genetic PC correction of the standard logistic regression GWAS analysis on all cause heart failure that we compared to our own method. 

## Directory: step7_adjust_HF_for_covariates_NN

 - `step7a_create_AE_phenotypes.py`: DO NOT RUN DIRECTLY. given a number of nodes in the autoencoder's first layer, the number of folds for cross validation, and the dropout probability, computes the average test accuracy of a trained autoencoder. The latent dimensions of this autoencoder comprise the autoencoder's (sometimes referred to as "NN" in path names) latent phenotypes. All 311 ICD codes and all cause heart failure were included as input.

 - `step7a_create_AE_phenotypes.sh`: Runs step7a_create_AE_phenotypes.py for 1000 nodes in the first layer, 5 fold cross validation, and dropout probabilities 0.05, ..., 0.5. From these, dropout probability 0.3 was selected for the next step. 

 - `step7b_create_best_phenotypes_normal_AE_0.3dop.py`: Uses the selected architecture to train a neural network without cross validation for the final model that produces autoencoder latent phenotypes. Model was trained twice prior to using it to produce the latent phenotypes. Training it the second time used resulting weights from the first time as starting points. 

 - `step7b_create_best_phenotypes_normal_AE_0.3dop.sh`: Run this to run step7b_create_best_phenotypes_normal_AE_0.3dop.py for training and retraining. 

 - `step7c_impute_missing_values.py`: DO NOT RUN DIRECTLY. This file applies MICE imputation to environmental factors' missing values. Features are filtered based on their strongest correlation with the environmental factors that missing values are being imputed for. A threshold for this correlation is specified via argparse parameter "nn". Also simulates different types of missingness to measure the imputation accuracy. 

 - `step7c_impute_missing_values.sh`: This file runs step7c_impute_missing_values.py to create latent features for nn = 0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, and 0.1. We ended up using 0.05 in the paper. 

 - `step7d_get_imputed_values_and_trasformed_variables.py`: essentially runs like step7c_impute_missing_values.py. The argparse nn parameter is not included, but instead set at 0.05. Also missingness is not simulated this time. Rather, missing values are imputed, and environmental factors with imputed missing values are produced. 

 - `step7d_get_imputed_values_and_trasformed_variables.sh`: Runs step7d_get_imputed_values_and_trasformed_variables.py. You could just run step7d_get_imputed_values_and_trasformed_variables.py directly, but it's runtime is such that submitting a job is less prone to disruption. 

 - `step7e_compute_network_shapley_values.py`: computes Shapley value contributions of ICD10 codes and all cause heart failure to a specifiedd autoencoder-based latent phenotype. 

 - `step7e_compute_network_shapley_values.sh`: applies step7e_compute_network_shapley_values.py to each such latent phenotype. 

 - `step7f_analyze_shapley_values.py`: generates dendrograms for correlations between shapley values with either high mean absolute values or high correlations to other features' shapley values. 

## Directory: step7_adjust_HF_for_covariates_PCA

 - `step0_compute_SV_p_values.py`: Analyzes which ICD10 codes' shapley values are most responsible for observed correlations between SNPs and latent phenotypes. Used by step7g_sub_phenotype_analysis.sh

 - `step7a_get_HF_ICD_codes_unrelated.py`: Imports phenotypes and relevant UK biobank fields corresponding to the eids of unrelated individuals determined in step 4.

 - `step7b_PCA_transform_the_data.py`: This file applies  PCA dimensionality reduction to 311 ICD10 codes and all cause heart failure. Also generates the 16th latent phenotype defined by the all-cause heart failure residuals. Adjusts all latent phenotypes by the principal components from step 6.

 - `step7c_impute_missing_values.py`: This file applies MICE imputation to environmental factors' missing values. Features are filtered based on their strongest correlation with the environmental factors that missing values are being imputed for. A threshold for this correlation is specified via argparse parameter "nn". Also simulates different types of missingness to measure the imputation accuracy. Rather than test multiple values, we selected nn = 0.05 for consistency and confirmed that the imputation accuracy is higher than mean imputation. 

 - `step7d_get_imputed_values_and_trasformed_variables.py`: essentially runs like step7c_impute_missing_values.py. The argparse nn parameter is not included, but instead set at 0.05. Also missingness is not simulated this time. Rather, missing values are imputed, and environmental factors with imputed missing values are produced. 

 - `step7e_compute_network_shapley_values.py`: computes Shapley value contributions of ICD10 codes and all cause heart failure to each PCA-based latent phenotype. 

 - `step7f_analyze_shapley_values.py`: generates dendrograms for correlations between shapley values with either high mean absolute values or high correlations to other features' shapley values. 

 - `step7g_sub_phenotype_analysis.sh`: Runs step0_compute_SV_p_values.py with different numbers of top ICD10 codes' shapley values to include in the subset analysis. 

## Directory: step8_get_imputed_ukb_samples

 - `step8.1_get_imputed_data_setup.py`: for each chromosome, creates a shell script for importing imputed SNPs on which GWAS is conducted.   
 
 - `step8.2_get_rsID_positions.py`: For each chromosome, creates a list of SNPs to removes for having a MAF that is too low or from containing too little information. 

 - `step8.3_get_eids.py`: creates a file specifying which unrelated individuals to keep.

 - `step8.4_get_imputed_data.sh`: Runs the bash scripts created in step8.1_get_imputed_data_setup.py

 - `step8.5_make_bims_tab_spaced.py`: does what the file name suggests. 

## Directory: step9_regress_phenotypes_against_SNPs_logistic_PCA

 - `step0_binary_HF_QTL_getter.py`: Regresses a subset of all SNPs against all cause heart failure with logistic regression and corrects with with genetic principal components. 

 - `step0_QQ_plot_getter.py`: Generates a QQ plot for all SNPs' real vs expected p values with respect to one latent phenotype.  

 - `step0_significant_QTL_getter.py`: Regresses all SNPs against one genetic PC corrected PCA based latent phenotype. Uses standard linear regression, EDGE, and target encoding. Must manually modify line 87 of the hard-code "env_name" to select which environmental factor is used in combination with target encoding to search for GxE effects.  

 - `step9a_get_genotype_metadata.sh`: gets SNPs' minor allele frequencies. Having these pre-computed makes the GWAS faster.

 - `step9b_get_significant_QTLs_setup.py`: creates one bash file per (latent phenotype, chromosome) pair for step0_binary_HF_QTL_getter.py to act on. Remember that the environmental factor must be hardcoded in step0_binary_HF_QTL_getter.py in replacement of line 87.

 - `step9c_get_significant_QTLs.sh`: runs all bash scripts produced by step9b_get_significant_QTLs_setup.py

 - `step9d_get_QQ_plots.sh`: applies step0_QQ_plot_getter.py to all latent phenotypes

 - `step9e_get_binary_HF_QTLs.sh`: creates bash scripts to run step0_binary_HF_QTL_getter.py on numerous smaller subsets of SNPs to parallelize the GWAS. 

 - `step9f_get_QQ_plots_normal_GWAS.sh`: applies step0_QQ_plot_getter.py to all cause heart failure. Some of the code also applying it to different binary phenotypes is obsolete.

## Directory: step9_regress_phenotypes_against_SNPs_NN

 - `step0_QQ_plot_getter.py`: Generates a QQ plot for all SNPs' real vs expected p values with respect to one latent phenotype.  

 - `step0_significant_QTL_getter.py`: Regresses all SNPs against one genetic PC corrected PCA based latent phenotype. Uses standard linear regression, EDGE, and target encoding. Must manually modify line 87 of the hard-code "env_name" to select which environmental factor is used in combination with target encoding to search for GxE effects.  

 - `step9a_get_genotype_metadata.sh`: gets SNPs' minor allele frequencies. Having these pre-computed makes the GWAS faster.

 - `step9b_get_significant_QTLs_setup.py`: creates one bash file per (latent phenotype, chromosome) pair for step0_binary_HF_QTL_getter.py to act on. Remember that the environmental factor must be hardcoded in step0_binary_HF_QTL_getter.py in replacement of line 87.

 - `step9c_get_significant_QTLs.sh`: runs all bash scripts produced by step9b_get_significant_QTLs_setup.py

 - `step9d_get_QQ_plots.sh`: applies step0_QQ_plot_getter.py to all latent phenotypes

## Directory: step9_regress_phenotypes_against_SNPs_PCA

 - `step0_QQ_plot_getter.py`: Generates a QQ plot for all SNPs' real vs expected p values with respect to one latent phenotype.  

 - `step0_significant_QTL_getter.py`: Regresses all SNPs against one genetic PC corrected PCA based latent phenotype. Uses standard linear regression, EDGE, and target encoding. Must manually modify line 87 of the hard-code "env_name" to select which environmental factor is used in combination with target encoding to search for GxE effects.  

 - `step9a_get_genotype_metadata.sh`: gets SNPs' minor allele frequencies. Having these pre-computed makes the GWAS faster. 

 - `step9b_get_significant_QTLs_setup.py`: creates one bash file per (latent phenotype, chromosome) pair for step0_binary_HF_QTL_getter.py to act on. Remember that the environmental factor must be hardcoded in step0_binary_HF_QTL_getter.py in replacement of line 87.

 - `step9c_get_significant_QTLs.sh`: runs all bash scripts produced by step9b_get_significant_QTLs_setup.py

 - `step9d_get_QQ_plots.sh`: applies step0_QQ_plot_getter.py to all latent phenotypes

## Directory: step10_get_significant_SNPs_logistic_PCA

 - `step0_compute_GxE_p_values.py`: given an rsID, chromosome, latent phenotype, and environmental factor as input, computes a permutation test p value for the pure GxE effect. 

 - `step0_filter_significant_SNPs_and_get_GxE_effects.py`: for a specified chromosome and environmental factor, for each latent phenotype, divides all SNP hits with those attributes into computationally tractable intervals. For each interval, selects independent SNPs with forward selection stepwise regression. Outputs files with independent SNP hits into the appropriate hits_QTL_[ENVIRONMENTAL FACTOR] directory. Also, for each independent SNP hit, creates a bash file in the "hits_GxE_p_vals_getters_[ENVIRONMENTAL FACTOR] directory that runs step0_compute_GxE_p_values.py

 - `step0_filter_significant_SNPs_and_get_GxE_effects_LAST_PART_ONLY.py`: If outputs files with independent SNP hits exist in the appropriate hits_QTL_[ENVIRONMENTAL FACTOR] directory, then for each independent SNP hit, creates a bash file in the "hits_GxE_p_vals_getters_[ENVIRONMENTAL FACTOR] directory that runs step0_compute_GxE_p_values.py

 - `step10a_get_significant_rsIDs.py`: Retreives all rsIDS corresponding to SNP hits with a nominal TRACE p value < 5E-8/16 (a bonferroni correction for the number of latent phenotypes)

 - `step10b_get_significant_SNPs.sh`: Retreives plink files for all SNPs corresponding to the rsIDs from the previous step

 - `step10c_filter_significant_SNPs_and_get_GxE_effects.sh`: applies step0_filter_significant_SNPs_and_get_GxE_effects.py to all chromosomes and environmental factors

 - `step10c_filter_significant_SNPs_and_get_GxE_effects_LAST_PART_ONLY.sh`: applies step0_filter_significant_SNPs_and_get_GxE_effects_LAST_PART_ONLY.py to all chromosomes and environmental factors 

 - `step10d_get_significant_GxE_p_values.sh`: applies step0_compute_GxE_p_values.py to all independently nominally significant rsIDs.  

 - `step10e_access_common_SNPs.py`: for the appropriate latent phenotype set (either PCA, logistic PCA, or the autoencoder), produces a list of main effects, gene by smoking interaction effects, gene by alcohol interaction effects, gene by gender interaction effects, and gene by exercise interaction effects. Also produces a list counting the number of each such effect for each latent phenotype. Also produces data (step10e_p_val_analysis.txt) for the p-value analysis in table 1b of the manuscript. Some code for machine learning exists beyond that point, but it is largely replaced with step10g.

 - `step10f_get_CV_folds.py`: Produces 30 index sets for outer training and outer validation folds that step10g uses to perform nested cross validation for the machine learning analysis. 

 - `step10g_get_CV_testing_accuracy.py`: for a specified set of the 30 index sets previously generated, computes 10 fold cross validation on the outer training set. Outputs optimal parameters for the gradient boosting model as well as accuracy on the outer validation set. 

 - `step10g_get_CV_testing_accuracy.sh`: applies step10g_get_CV_testing_accuracy.py to all 30 folds. 

## Directory: step10_get_significant_SNPs_NN

 - `step0_compute_GxE_p_values.py`: given an rsID, chromosome, latent phenotype, and environmental factor as input, computes a permutation test p value for the pure GxE effect. 

 - `step0_filter_significant_SNPs_and_get_GxE_effects.py`: for a specified chromosome and environmental factor, for each latent phenotype, divides all SNP hits with those attributes into computationally tractable intervals. For each interval, selects independent SNPs with forward selection stepwise regression. Outputs files with independent SNP hits into the appropriate hits_QTL_[ENVIRONMENTAL FACTOR] directory. Also, for each independent SNP hit, creates a bash file in the "hits_GxE_p_vals_getters_[ENVIRONMENTAL FACTOR] directory that runs step0_compute_GxE_p_values.py

 - `step0_filter_significant_SNPs_and_get_GxE_effects_LAST_PART_ONLY.py`: If outputs files with independent SNP hits exist in the appropriate hits_QTL_[ENVIRONMENTAL FACTOR] directory, then for each independent SNP hit, creates a bash file in the "hits_GxE_p_vals_getters_[ENVIRONMENTAL FACTOR] directory that runs step0_compute_GxE_p_values.py

 - `step10a_get_significant_rsIDs.py`: Retreives all rsIDS corresponding to SNP hits with a nominal TRACE p value < 5E-8/16 (a bonferroni correction for the number of latent phenotypes)

 - `step10b_get_significant_SNPs.sh`: Retreives plink files for all SNPs corresponding to the rsIDs from the previous step

 - `step10c_filter_significant_SNPs_and_get_GxE_effects.sh`: applies step0_filter_significant_SNPs_and_get_GxE_effects.py to all chromosomes and environmental factors

 - `step10c_filter_significant_SNPs_and_get_GxE_effects_LAST_PART_ONLY.sh`: applies step0_filter_significant_SNPs_and_get_GxE_effects_LAST_PART_ONLY.py to all chromosomes and environmental factors  

 - `step10d_get_significant_GxE_p_values.sh`: applies step0_compute_GxE_p_values.py to all independently nominally significant rsIDs.  

 - `step10e_access_common_SNPs.py`: for the appropriate latent phenotype set (either PCA, logistic PCA, or the autoencoder), produces a list of main effects, gene by smoking interaction effects, gene by alcohol interaction effects, gene by gender interaction effects, and gene by exercise interaction effects. Also produces a list counting the number of each such effect for each latent phenotype. Also produces data (step10e_p_val_analysis.txt) for the p-value analysis in table 1b of the manuscript. Some code for machine learning exists beyond that point, but it is largely replaced with step10g. 

 - `step10f_get_CV_folds.py`: Produces 30 index sets for outer training and outer validation folds that step10g uses to perform nested cross validation for the machine learning analysis. 

 - `step10g_get_CV_testing_accuracy.py`: for a specified set of the 30 index sets previously generated, computes 10 fold cross validation on the outer training set. Outputs optimal parameters for the gradient boosting model as well as accuracy on the outer validation set.  

 - `step10g_get_CV_testing_accuracy.sh`: applies step10g_get_CV_testing_accuracy.py to all 30 folds. 

## Directory: step10_get_significant_SNPs_PCA

 - `step0_compute_GxE_p_values.py`: given an rsID, chromosome, latent phenotype, and environmental factor as input, computes a permutation test p value for the pure GxE effect. 

 - `step0_filter_significant_SNPs_and_get_GxE_effects.py`: for a specified chromosome and environmental factor, for each latent phenotype, divides all SNP hits with those attributes into computationally tractable intervals. For each interval, selects independent SNPs with forward selection stepwise regression. Outputs files with independent SNP hits into the appropriate hits_QTL_[ENVIRONMENTAL FACTOR] directory. Also, for each independent SNP hit, creates a bash file in the "hits_GxE_p_vals_getters_[ENVIRONMENTAL FACTOR] directory that runs step0_compute_GxE_p_values.py

 - `step0_filter_significant_SNPs_and_get_GxE_effects_LAST_PART_ONLY.py`: If outputs files with independent SNP hits exist in the appropriate hits_QTL_[ENVIRONMENTAL FACTOR] directory, then for each independent SNP hit, creates a bash file in the "hits_GxE_p_vals_getters_[ENVIRONMENTAL FACTOR] directory that runs step0_compute_GxE_p_values.py

 - `step10a_get_significant_rsIDs.py`: Retreives all rsIDS corresponding to SNP hits with a nominal TRACE p value < 5E-8/16 (a bonferroni correction for the number of latent phenotypes). 

 - `step10b_get_significant_SNPs.sh`: Retreives plink files for all SNPs corresponding to the rsIDs from the previous step

 - `step10c_filter_significant_SNPs_and_get_GxE_effects.sh`: applies step0_filter_significant_SNPs_and_get_GxE_effects.py to all chromosomes and environmental factors

 - `step10c_filter_significant_SNPs_and_get_GxE_effects_LAST_PART_ONLY.sh`: applies step0_filter_significant_SNPs_and_get_GxE_effects_LAST_PART_ONLY.py to all chromosomes and environmental factors  

 - `step10d_get_significant_GxE_p_values.sh`: applies step0_compute_GxE_p_values.py to all independently nominally significant rsIDs.  

 - `step10e_access_common_SNPs.py`: for the appropriate latent phenotype set (either PCA, logistic PCA, or the autoencoder), produces a list of main effects, gene by smoking interaction effects, gene by alcohol interaction effects, gene by gender interaction effects, and gene by exercise interaction effects. Also produces a list counting the number of each such effect for each latent phenotype. Also produces data (step10e_p_val_analysis.txt) for the p-value analysis in table 1b of the manuscript. Some code for machine learning exists beyond that point, but it is largely replaced with step10g.

 - `step10f_get_CV_folds.py`: Produces 30 index sets for outer training and outer validation folds that step10g uses to perform nested cross validation for the machine learning analysis. 

 - `step10g_get_CV_testing_accuracy.py`: for a specified set of the 30 index sets previously generated, computes 10 fold cross validation on the outer training set. Outputs optimal parameters for the gradient boosting model as well as accuracy on the outer validation set. 

 - `step10g_get_CV_testing_accuracy.sh`: applies step10g_get_CV_testing_accuracy.py to all 30 folds. 

## Directory: step11_analyze_complete_dataset

 - `step11a_merge_rsID_output.py`: generates table 1a from the manuscript. Also generates five files of the form rsIDs_*_effects.txt, and five of the form rsIDs_*_effects_pvals.txt. The astrisk is a variable part of the filename specifying the effect type.

 - `step11a_merge_rsID_output.py IMPORTANT OUTPUT NOTE`: The (rsIDs_*_effects.txt, rsIDs_*_effects_pvals.txt) file pair for each respective effect type is used as input for FUMA's snp2gene tool. rsIDs_*_effects_pvals.txt is input for "GWAS summary statistics", and rsIDs_*_effects.txt is input for Pre-defined lead SNPs. We unchecked the "Identify additional independent lead SNPs" box, input 380000 for a total sample size (not relevant because we're not using this for statistical analysis), and set the maximum distance to 500KB in the positional mapping box. Genes farther than 300KB from SNP hits were filtered out later. We have renamed the relevant FUMA output files to annov_alcohol.txt, annov_gender.txt, annov_smoking.txt, and annov_main.txt. 

 - `step11b_get_enrichment_all.R`: DO NOT RUN THE FILE ALL AT ONCE. There are clearly outlined manual steps that must be taken to run the miRNA enrichment analysis. With that said, this script uses renamed FUMA output files beginning with "annov" as input. Also uses rsIDs_*_effects.txt files as input. Generates table 2 and table S6. The end result of the manual part is to produce the (manually renamed) file miEAA.csv, which lists kegg pathways enriched for miRNAs that were previously enriched for our SNP hits. 

 - `step11c_get_chr_seperated_lists.py`: For each SNP hit, makes a single file that step11d_ldlink.sh remotely sends as input to https://ldlink.nih.gov/?tab=ldtrait

 - `step11d_ldlink.sh`: sends each output from step11c_get_chr_seperated_lists.py as input to https://ldlink.nih.gov/?tab=ldtrait. For each SNP, finds all prior GWAS SNP hits from other research that were in LD with our SNP and downloads relevant data about those SNPs. Thus, each input file of a single SNP either corresponds to nothing at all or an output file with one or more prior GWAS SNP hits in LD with our input SNP hit.  

 - `step11e_make_SNP_tables_LD_pruning.py`: Uses the output from step11d_ldlink.sh as input to make table 1b. It counts the number of independent SNP hits (defined as having an R squared LD of less than 0.8) that are prior GWAS hits for phenotypes that are relevant to all cause heart failure. Our list of search terms that we believe are relevant to all cause heart failure was manually created and named "possible_AHF_terms" in the file. They are substrings of names for biological topics that are closely related to cardiovascular dysfunction. 

 - `step11f_miRNA_enrichment_analysis.R`: DO NOT RUN THE FILE ALL AT ONCE. There are clearly outlined manual steps. This file computes p-values for miRNA associated SNP hits being enriched for genic SNPs. Also generates table 3 and table S7. The content in genevestigator_hits was created by using genevestigator's interface, and the methods are detailed in "genevestigator_hits_methods". This comprises the gene-smoking associations used in tables 3 and S7. There is also a manual step used to get the gene-disease associations.

 - `step11g_make_heritability_figure.py`: uses output from the step10g_get_CV_testing_accuracy.sh files in all three directories (PCA, logistic_PCA, NN) to make table S5b (best gradient boosting classification inner fold cross validation hyperparameters) and figure 3a. 

 - `step11h_confirm_model_consistency.py`: generates table 1c and table 1d as described in the manuscript's methods. 

 - `step11i_get_independent_normal_GWAS_SNPs.py`: finds any SNP hits for standard logistic regression GWAS that are independently significant from rs1906609, rs7857118, rs12627426, rs73839819, rs2234962, and rs12138073. Since these are the SNPs analyzed by the paper defining all cause heart failure, the code assumes these are the top hits by default.

 - `step11j_make_figure3b.sh`: simply makes the files necessary for producing figure 3b. SNPs included in the analysis include the six default SNP and one other independently significant SNP that the previous step had found (rs73188900). 

 - `step11k_make_table1.R`: Takes the output from the previous step to make figure 3b. Also makes table 1. 

 - `step11l_make_supp_figs.py`: Makes tables S1, S2, S3, and S4. Makes figures S1 and S2. Makes some input for figure S3a. 

 - `step11m_finish_fig_S3a.R`: Finishes figure S3a. 


