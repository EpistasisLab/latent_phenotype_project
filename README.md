# latent_phenotype_project 

NOTE: File names are listed and described in the order that they are supposed to be run. 

NOTE: Any file name containing "step0" is not meant to be run directly. Rather, some other file in the sequence calls that file to run when necessary.

NOTE: Descriptions for files that are called by other files but do not have "step0" in the name start with "DO NOT RUN DIRECTLY"

NOTE: Directories are ordered from top to bottom as the sequence in which they should be run. 

## Directory: step1_get_phenotypes_complex

 - `step1a_get_UKB_phenotypes.py`: Imports and processes UK biobank (UKB) data fields. Handles multiple measurements, modifies smoking pack years, calculates annual alcohol consumption, and binarizes categorical variables.

 - `step1a_library.py`: Functions used by `step1a_get_UKB_phenotypes.py`.

## Directory: step1_get_phenotypes_simple

 - `step1a_import_HF_ICD_codes_and_covariates.py`: Imports additional UKB data fields related to heart failure; binarizes and renames as needed.

 - `step1b_process_HF_ICD_codes_and_covariates.py`: Binarizes ICD codes, creates the all cause heart failure binary phenotype ('AHF' in text, 'any_HF' in code), and removes ICD codes not significantly correlated to AHF.

## Directory: step2_get_UKB_samples

 - `step2a_make_UKB_sample_getters.py`:  Generates bash scripts for importing individual' non-imputed genotype data.

 - `step2b_get_UKB_samples.sh`: Executes bash scripts from the previous step to collect genotype samples.

## Directory: step3_merge_chr_and_remove_quitters

 - `step3a_get_people_who_quit.py`: generates a list of individuals who opted out of UKB data use (from other manually typed lists).

 - `step3b_merge_datasets.sh`: Merges plink files from Step 2; removes opt-out individuals.

 - `step3c_find_SNP_cutoffs.py`: Plots SNPs' minor allele frequencies, missingness, and HWE p-values (exploratory only).

 - `step3d_remove_bad_SNPs.sh`: Removes SNPs exceeding thresholds for aforementioned quantities.

 - `step3e_find_sample_cutoffs.py`: Plots individuals' heterozygosity, average SNP missingness, and X chromosome heterozygosity scores (exploratory only).

 - `step3f_remove_bad_samples.sh`: Removes individuals exceeding thresholds for aforementioned quantities.

## Directory: step4_remove_relatives

 - `step4a_divide_eids_into_subsets_prep.py`: Splits remaining eids from Step 3 into 10 equal subsets. 

 - `step4b_divide_eids_into_subsets.sh`: Partitions plink files based on the 10 eid subsets; removes remainder eids.

 - `step4c_make_king_subjobs.py`: Generates a bash script that runs KING for each SNP subset and each pair of SNP subsets.

 - `step4d_run_king_subjobs.sh`: Runs KING through the previously generated bash scripts. King outputs all pairs of related individuals with 3rd degree relatedness or more. 

 - `step4e_get_unrelated_eids.py`: Prunes dataset to keep only 4th degree or less related individuals. The number of all-cause heart failure cases is prioritized over the total sample size. Both are maximized with that in mind. 

 - `step4f_remove_relatives.sh`: Excludes relatives based on prior steps. From these unrelated individuals, creates a pruned SNP dataset such that all SNPs are in low LD for genetic principle component computation. Also keeps the non-pruned dataset.

## Directory: step5_verify_sample_unrelatedness

- `step5a_divide_eids_into_subsets_prep.py`: Splits eids from Step 3 into 10 equal subsets. 

- `step5b_divide_eids_into_subsets.sh`: Partitions plink files into 10 subsets corresponding to the eid subsets.

- `step5c_make_king_subjobs.py`: Generates a bash script to run KING for each SNP subset and each pair of SNP subsets.

- `step5d_run_king_subjobs.sh`: Executes the KING subjob scripts.

- `step5e_confirm_unrelatedness.py`: Confirms that no related individuals remain post-Step 4. This validates that KING was used correctly.

## Directory: step6_PCA

- `check_LD.sh`: Computes LD between SNP pairs in pruned SNP set from Step 4.

- `check_LD.py`: Validates that plink was used correctly, so no SNP pairs exceed the LD R^2 threshold from Step 4.
 
- `step0_run_PCA.par`: Specifies parameters for Eigensoft's PCA computation.

- `step6a_make_PCA_documents.py`: Prepares input files for Eigensoft's PCA.

- `step6b_run_PCA.sh`: Executes PCA computation using Eigensoft.

## Directory: step7_adjust_HF_for_covariates_logistic_PCA

 - `step7a_get_HF_ICD_codes_unrelated.py`: Imports phenotypes and relevant UK Biobank fields for unrelated individuals from Step 4.

 - `step7b_logPCA_transform_the_data.py`: DO NOT RUN DIRECTLY. Applies logistic PCA to 311 ICD10 codes and heart failure. Creates latent phenotypes (k specified via argparse) and adjusts them using PCs from Step 6.

 - `step7b_logPCA_transform_the_data.sh`: Executes the above for k = 1 to 20; paper uses k=15.

 - `step7c_impute_missing_values.py`: DO NOT RUN DIRECTLY. Applies MICE imputation to environmental factors with missing values. Correlation threshold for feature selection specified via "nn".

 - `step7c_impute_missing_values.sh`: Runs the above for various nn values; paper uses nn=0.05. 

 - `step7d_get_imputed_values_and_trasformed_variables.py`: Similar to above, but sets nn at 0.05 and does not simulate missingness. Outputs imputed environmental factors.

 - `step7d_get_imputed_values_and_trasformed_variables.sh`: Executes the above, recommended for job submission due to long runtime.

 - `step7e_get_PCs_effs.py`: Outputs logistic regression beta coefficients for genetic PCs vs AHF, used to correct the standard logistic regression GWAS against AHF. 

## Directory: step7_adjust_HF_for_covariates_NN

- `step7a_create_AE_phenotypes.py`: DO NOT RUN DIRECTLY. Computes autoencoder test accuracy based on layer nodes, cross-validation folds, and dropout rate. Generates latent phenotypes.

- `step7a_create_AE_phenotypes.sh`: Runs the above for specified nodes, folds, and dropout rates; paper uses 0.3 dropout.

- `step7b_create_best_phenotypes_normal_AE_0.3dop.py`: Trains final autoencoder model twice, using the first run's weights as a starting point for the second.

- `step7b_create_best_phenotypes_normal_AE_0.3dop.sh`: Executes the above model training.

- `step7c_impute_missing_values.py`: DO NOT RUN DIRECTLY. Uses MICE to impute missing environmental data, selecting features based on correlation "nn".

- `step7c_impute_missing_values.sh`: Executes imputation for various "nn"; paper uses nn=0.05.

- `step7d_get_imputed_values_and_transformed_variables.py`: Similar to above, fixes "nn" at 0.05 and outputs imputed factors.

- `step7d_get_imputed_values_and_transformed_variables.sh`: Executes the above, recommended for job submission.

- `step7e_compute_network_shapley_values.py`: Calculates Shapley values for each latent phenotype.

- `step7e_compute_network_shapley_values.sh`: Applies the above to all latent phenotypes.

- `step7f_analyze_shapley_values.py`: Produces dendrograms based on correlations of high-impact or high-correlation Shapley values.

## Directory: step7_adjust_HF_for_covariates_PCA

- `step0_compute_SV_p_values.py`: Identifies key ICD10 codes affecting SNP-latent phenotype correlations. Used by `step7g_sub_phenotype_analysis.sh`.

- `step7a_get_HF_ICD_codes_unrelated.py`: Fetches unrelated individuals' phenotypes and UK Biobank fields from Step 4.

- `step7b_PCA_transform_the_data.py`: Performs PCA on 311 ICD10 codes and all-cause heart failure, adjusting latent phenotypes with PCs from Step 6.

- `step7c_impute_missing_values.py`: Applies MICE imputation to environmental factors, using a fixed "nn" of 0.05. Confirmed to outperform mean imputation.

- `step7d_get_imputed_values_and_transformed_variables.py`: Similar to `step7c`, but with "nn" fixed at 0.05 and no missingness simulation. Outputs imputed factors.

- `step7e_compute_network_shapley_values.py`: Calculates Shapley values for ICD10 codes and heart failure in PCA-based latent phenotypes.

- `step7f_analyze_shapley_values.py`: Generates dendrograms based on Shapley value correlations.

- `step7g_sub_phenotype_analysis.sh`: Executes `step0_compute_SV_p_values.py` for varying numbers of top-contributing ICD10 codes.

## Directory: step8_get_imputed_ukb_samples

- `step8.1_get_imputed_data_setup.py`: Generates shell scripts for importing GWAS-targeted, imputed SNPs per chromosome.

- `step8.2_get_rsID_positions.py`: Lists SNPs to remove per chromosome based on low MAF or insufficient information.

- `step8.3_get_eids.py`: copies a file identifying which unrelated individuals to retain in the analysis.

- `step8.4_get_imputed_data.sh`: Executes shell scripts from `step8.1`.

- `step8.5_make_bims_tab_spaced.py`: Converts BIM files to tab-spaced format.

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

