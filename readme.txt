# latent_phenotype_project 

NOTE: File names are listed and described in the order that they are supposed to be run. 

NOTE: Any file name containing "step0" is not meant to be run directly. Rather, some other file in the sequence calls that file to run when necessary.

NOTE: Descriptions for files that are called by other files but do not have "step0" in the name start with "DO NOT RUN DIRECTLY"

NOTE: Descriptions starting with "IMPORTANT; run in sections, not all at once" have sections that need to be completed manually in the code. You may run the code in subsections that do not overlap the sections to be completed manually. 

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

 - `step0_binary_HF_QTL_getter.py`: Performs logistic regression on a subset of SNPs against all-cause heart failure, corrected by genetic PCs. 

 - `step0_QQ_plot_getter.py`: Generates a QQ plot for all SNPs' real vs expected p values with respect to one latent phenotype.  

 - `step0_significant_QTL_getter.py`: Regresses SNPs against a genetic PC-corrected latent phenotype using linear regression, EDGE, and target encoding. Manually set "env_name" on line 87 for GxE effects. Examined Env_name values include 'pack-years', 'annual-consumption', ['874-average', '894-average', '914-average'] (a python list), and '22001-0.0'

 - `step9a_get_genotype_metadata.sh`: gets SNPs' minor allele frequencies. Having these pre-computed makes the GWAS faster.

 - `step9b_get_significant_QTLs_setup.py`: Creates bash scripts for each (latent phenotype, chromosome) pair. Remember to update line 87 in `step0_binary_HF_QTL_getter.py` to change the environmental factor.

 - `step9c_get_significant_QTLs.sh`: Executes scripts generated by the previous step.

 - `step9d_get_QQ_plots.sh`: Applies `step0_QQ_plot_getter.py` to all latent phenotypes.

 - `step9e_get_binary_HF_QTLs.sh`: creates bash scripts to run `step0_binary_HF_QTL_getter.py` on all SNP subsets in parallel. 

 - `step9f_get_QQ_plots_normal_GWAS.sh`: applies `step0_QQ_plot_getter.py` to all cause heart failure.

## Directory: step9_regress_phenotypes_against_SNPs_NN

 - `step0_QQ_plot_getter.py`: Generates a QQ plot for all SNPs' real vs expected p values with respect to one latent phenotype.  

 - `step0_significant_QTL_getter.py`: Regresses SNPs against a genetic PC-corrected latent phenotype using linear regression, EDGE, and target encoding. Manually set "env_name" on line 87 for GxE effects. Examined Env_name values include 'pack-years', 'annual-consumption', ['874-average', '894-average', '914-average'] (a python list), and '22001-0.0'

 - `step9a_get_genotype_metadata.sh`: gets SNPs' minor allele frequencies. Having these pre-computed makes the GWAS faster.

 - `step9b_get_significant_QTLs_setup.py`: Creates bash scripts for each (latent phenotype, chromosome) pair. Remember to update line 87 in `step0_binary_HF_QTL_getter.py` to change the environmental factor.

 - `step9c_get_significant_QTLs.sh`: Executes scripts generated by the previous step.

 - `step9d_get_QQ_plots.sh`: Applies `step0_QQ_plot_getter.py` to all latent phenotypes.

## Directory: step9_regress_phenotypes_against_SNPs_PCA

 - `step0_QQ_plot_getter.py`: Generates a QQ plot for all SNPs' real vs expected p values with respect to one latent phenotype.  

 - `step0_significant_QTL_getter.py`: Regresses SNPs against a genetic PC-corrected latent phenotype using linear regression, EDGE, and target encoding. Manually set "env_name" on line 87 for GxE effects. Examined Env_name values include 'pack-years', 'annual-consumption', ['874-average', '894-average', '914-average'] (a python list), and '22001-0.0'

 - `step9a_get_genotype_metadata.sh`: gets SNPs' minor allele frequencies. Having these pre-computed makes the GWAS faster. 

 - `step9b_get_significant_QTLs_setup.py`: Creates bash scripts for each (latent phenotype, chromosome) pair. Remember to update line 87 in `step0_binary_HF_QTL_getter.py` to change the environmental factor.

 - `step9c_get_significant_QTLs.sh`: Executes scripts generated by the previous step.

 - `step9d_get_QQ_plots.sh`: Applies `step0_QQ_plot_getter.py` to all latent phenotypes.

## Directory: step10_get_significant_SNPs_logistic_PCA

 - `step0_compute_GxE_p_values.py`: given an rsID, chromosome, latent phenotype, and environmental factor as input, computes a permutation test p value for the pure GxE effect. 

 - `step0_filter_significant_SNPs_and_get_GxE_effects.py`: For each chromosome and environmental factor, for each latent phenotype, segments SNP hits into intervals, and selects independently nominally significant SNPs. For each independent SNP hit, prepares a bash file to run `step0_compute_GxE_p_values.py`

 - `step0_filter_significant_SNPs_and_get_GxE_effects_LAST_PART_ONLY.py`: If independent SNP hits files are available, prepares bash files to run `step0_compute_GxE_p_values.py as described previously.

 - `step10a_get_significant_rsIDs.py`: Retreives all rsIDS corresponding to SNP hits with a nominal TRACE p value < 5E-8/16 (16 is a bonferroni correction for the number of latent phenotypes)

 - `step10b_get_significant_SNPs.sh`: Retreives plink files for all SNPs from the previous step

 - `step10c_filter_significant_SNPs_and_get_GxE_effects.sh`: applies `step0_filter_significant_SNPs_and_get_GxE_effects.py` to all chromosomes and environmental factors

 - `step10c_filter_significant_SNPs_and_get_GxE_effects_LAST_PART_ONLY.sh`: applies `step0_filter_significant_SNPs_and_get_GxE_effects_LAST_PART_ONLY.py` to all chromosomes and environmental factors 

 - `step10d_get_significant_GxE_p_values.sh`: applies `step0_compute_GxE_p_values.py` to all independently nominally significant rsIDs.  

 - `step10e_access_common_SNPs.py`: Generates lists of main and interaction effects for various environmental factors with logistic PCA latent phenotypes, counts of each effect, and data for table 1a (step10e_p_val_analysis.txt). Machine learning is deferred to step10g. 

 - `step10f_get_CV_folds.py`: Generates 30 outer training/validation index sets for nested cross-validation.

 - `step10g_get_CV_testing_accuracy.py`: For 1 of the 30 outer index sets, conducts 10-fold cross-validation on the training set. Reports optimal model parameters and validation accuracy.

 - `step10g_get_CV_testing_accuracy.sh`: applies `step10g_get_CV_testing_accuracy.py` to all 30 outer training/validation index sets.

## Directory: step10_get_significant_SNPs_NN

 - `step0_compute_GxE_p_values.py`: given an rsID, chromosome, latent phenotype, and environmental factor as input, computes a permutation test p value for the pure GxE effect. 

 - `step0_filter_significant_SNPs_and_get_GxE_effects.py`: For each chromosome and environmental factor, for each latent phenotype, segments SNP hits into intervals, and selects independently nominally significant SNPs. For each independent SNP hit, prepares a bash file to run `step0_compute_GxE_p_values.py`

 - `step0_filter_significant_SNPs_and_get_GxE_effects_LAST_PART_ONLY.py`: If independent SNP hits files are available, prepares bash files to run `step0_compute_GxE_p_values.py as described previously.

 - `step10a_get_significant_rsIDs.py`: Retreives all rsIDS corresponding to SNP hits with a nominal TRACE p value < 5E-8/16 (16 is a bonferroni correction for the number of latent phenotypes)

 - `step10b_get_significant_SNPs.sh`: Retreives plink files for all SNPs from the previous step

 - `step10c_filter_significant_SNPs_and_get_GxE_effects.sh`: applies `step0_filter_significant_SNPs_and_get_GxE_effects.py` to all chromosomes and environmental factors

 - `step10c_filter_significant_SNPs_and_get_GxE_effects_LAST_PART_ONLY.sh`: applies `step0_filter_significant_SNPs_and_get_GxE_effects_LAST_PART_ONLY.py` to all chromosomes and environmental factors  

 - `step10d_get_significant_GxE_p_values.sh`: applies `step0_compute_GxE_p_values.py` to all independently nominally significant rsIDs.  

 - `step10e_access_common_SNPs.py`: Generates lists of main and interaction effects for various environmental factors with NN latent phenotypes, counts of each effect, and data for table 1a (step10e_p_val_analysis.txt). Machine learning is deferred to step10g. 

 - `step10f_get_CV_folds.py`: Generates 30 outer training/validation index sets for nested cross-validation.

 - `step10g_get_CV_testing_accuracy.py`: For 1 of the 30 outer index sets, conducts 10-fold cross-validation on the training set. Reports optimal model parameters and validation accuracy.

 - `step10g_get_CV_testing_accuracy.sh`: applies `step10g_get_CV_testing_accuracy.py` to all 30 outer training/validation index sets.

## Directory: step10_get_significant_SNPs_PCA

 - `step0_compute_GxE_p_values.py`: given an rsID, chromosome, latent phenotype, and environmental factor as input, computes a permutation test p value for the pure GxE effect. 

 - `step0_filter_significant_SNPs_and_get_GxE_effects.py`: For each chromosome and environmental factor, for each latent phenotype, segments SNP hits into intervals, and selects independently nominally significant SNPs. For each independent SNP hit, prepares a bash file to run `step0_compute_GxE_p_values.py`

 - `step0_filter_significant_SNPs_and_get_GxE_effects_LAST_PART_ONLY.py`: If independent SNP hits files are available, prepares bash files to run `step0_compute_GxE_p_values.py as described previously.

 - `step10a_get_significant_rsIDs.py`: Retreives all rsIDS corresponding to SNP hits with a nominal TRACE p value < 5E-8/16 (16 is a bonferroni correction for the number of latent phenotypes). 

 - `step10b_get_significant_SNPs.sh`: Retreives plink files for all SNPs from the previous step

 - `step10c_filter_significant_SNPs_and_get_GxE_effects.sh`: applies `step0_filter_significant_SNPs_and_get_GxE_effects.py` to all chromosomes and environmental factors

 - `step10c_filter_significant_SNPs_and_get_GxE_effects_LAST_PART_ONLY.sh`: applies `step0_filter_significant_SNPs_and_get_GxE_effects_LAST_PART_ONLY.py` to all chromosomes and environmental factors  

 - `step10d_get_significant_GxE_p_values.sh`: applies `step0_compute_GxE_p_values.py` to all independently nominally significant rsIDs.  

 - `step10e_access_common_SNPs.py`: Generates lists of main and interaction effects for various environmental factors with PCA latent phenotypes, counts of each effect, and data for table 1 (step10e_p_val_analysis.txt). Machine learning is deferred to step10g.

 - `step10f_get_CV_folds.py`: Generates 30 outer training/validation index sets for nested cross-validation.

 - `step10g_get_CV_testing_accuracy.py`: For 1 of the 30 outer index sets, conducts 10-fold cross-validation on the training set. Reports optimal model parameters and validation accuracy.

 - `step10g_get_CV_testing_accuracy.sh`: applies `step10g_get_CV_testing_accuracy.py` to all 30 outer training/validation index sets.

## Directory: step11_analyze_complete_dataset

### Code Files

 - `step11a_merge_rsID_output.py`: Generates Table 1b and effect-specific rsID file pairs (rsIDs_*_effects.txt and rsIDs_*_effects_pvals.txt, where * is the effect type)  for FUMA input. 

 - `Important Note on FUMA details`: For each rsID file pair, Inputs both files into FUMA with a placeholder total sample size of 380000. It is merely a required input that does not effect the distances between SNP hits and genes, thereby making it irrelevant. Initially uses a 500KB gene-SNP distance, later refined to 300KB for precision with small gene count reduction. Skips gene by exercise interactions due to limited hits. Renames FUMA output files as annov_*.txt. for clarity, where * = (main, smoking, gender, alcohol).

 - `step11b_get_enrichment_all.R`: IMPORTANT; run in sections, not all at once. Intermediate manual steps generate miEAA.csv, which lists KEGG-enriched miRNA pathways. Uses FUMA outputs and rsIDs_*_effects.txt as input. Generates Tables 2 and S6. 

 - `step11c_get_chr_seperated_lists.py`:  makes one file per SNP hit for input into the LDlink tool. 

 - `step11d_ldlink.sh`:  Sends files from previous step to LDlink, fetching data on prior GWAS hits in LD with each SNP hit. 

 - `step11e_make_SNP_tables_LD_pruning.py`: Creates Table 1a, counting independent GWAS SNP hits that output from step11d_ldlink.sh finds are related to AHF. The term list "possible_AHF_terms" was manually curated and validated through substring matching of terms related to cardiovascular dysfunction in the output from step11d_ldlink.sh. 

 - `step11f_miRNA_enrichment_analysis.R`: IMPORTANT; run in sections, not all at once. Computes p-values for enrichment of genic SNP hits in genes that enrich miRNA. Produces Tables 3 and S7. Manual steps are clearly outlined for gene-disease associations. Also uses data from Genevestigator ("see genevestigator_hits_methods folder").

 - `step11g_make_heritability_figure.py`: Generates Table S5b and Figure 3a using outputs from step10g_get_CV_testing_accuracy.sh across PCA, logistic_PCA, and NN directories.

 - `step11h_confirm_model_consistency.py`: Generates Tables 1c and 1d as per manuscript methods.
  
 - `step11i_get_independent_normal_GWAS_SNPs.py`: Identifies SNP hits from standard logistic GWAS that are independent of six default AHF-defining SNPs.

 - `step11j_make_figure3b.sh`: Prepares files for Figure 3b, including six default SNPs and one additional independent SNP (rs73188900) found in the previous step.

 - `step11k_make_table1.R`: Creates Figure 3b and Table 1 using the output from the previous step. Also generates figure 1.

 - `step11l_make_supp_figs.py`: Produces Tables S1-S4 and Figures S1, S2, and partial input for Figure S3a.

 - `step11m_finish_fig_S3a.R`: Completes Figure S3a.

### Intermediate Files

 - `annov_alcohol.txt`: relevant FUMA output when using `rsIDs_GxAlcohol_effects.txt` and `rsIDs_GxAlcohol_effects_pvals.txt` as input. Originally named "annov.txt".

 - `annov_gender.txt`: relevant FUMA output when using `rsIDs_GxGender_effects.txt` and `rsIDs_GxGender_effects_pvals.txt` as input. Originally named "annov.txt".

 - `annov_main.txt`: relevant FUMA output when using `rsIDs_main_effects.txt` and `rsIDs_main_effects_pvals.txt` as input. Originally named "annov.txt".

 - `annov_smoking.txt`: relevant FUMA output when using `rsIDs_GxSmoking_effects.txt` and `rsIDs_GxSmoking_effects_pvals.txt` as input. Originally named "annov.txt".

 - `disease_gene_associations_logistic_PCA_smoking.tsv`: Output from entering "CASZ1::AKR7A3::PPIE::HIVEP3::SSBP3::ST6GALNAC3::DPYD::NOS1AP::DDR2::SLC9A2::THSD7B::RAPGEF4::GPR155::ALS2::CHL1::IL5RA::THRB::CACNA2D3::CADPS::PRICKLE2::ROBO1::GAP43::LSAMP::AGTR1::MED12L::SERPINI1::MECOM::PARL::ST6GAL1::MUC4::EVC2::NPNT::HPGD::TENM3::LPCAT1::DROSHA::PRLR::GHR::EDIL3::CHSY3::FAM53C::EXOC2::ATXN1::BTBD9::TFEB::RUNX2::ESR1::TMEM242::PDE10A::COX19::ICA1::CREB5::PKD1L1::CALN1::CALCR::RINT1::ATXN7L1::DGKI::DPP6::TNFRSF10B::ADRA1A::ASPH::HNF4G::UBR5::TBC1D31::TMEM65::SH3GL2::FBP1::PALM2-AKAP2::GARNL3::FAM107B::CXCL12::MYPN::CDH23::CCDC147::SORCS3::RBM20::GRK5::ATE1::SYT9::SBF2::SLC22A8::TENM4::MAML2::NTM::TMTC1::ANO4::ATP12A::LRCH1::HTR2A::KLHL1::GPC6::FNTB::GALNT16::SMOC1::MAP3K9::KCNK13::C14orf159::ITPK1::UNC79::EML1::RYR3::WDR72::UNC13C::ARNT2::KIAA1199::PDE8A::SLCO3A1::CDIP1::NXN::PIK3R5::CA10::GAA::PIEZO2::BCL2::CLEC4M::ZNF627::SLC1A6::KDELR1::MACROD2::PCSK2::STX16-NPEPL1::CDH4::RIPK4::POF1B" into https://www.disgenet.org/search. Output file was originally named "54897__22977__10450__59269__23648__256435__1806__9722__4921__6549__80731__11069__151556__57679__10752__3568__7068__55799__8618__166336__6091__2596__4045__185__116931__5274__2122__55486__6480__4585__132884__2557.tsv" and was downloaded on 9/12/2023. 

 - `disease_gene_associations_NN_smoking.tsv`: Output from entering "EPHA8::PGM1::SCCPDH::CNTNAP5::PLEKHM3::CNTN3::MAP1B::GRK6::GLP1R::COL19A1::SCARA5::ANP32B::NAV2::KSR2::TMEM132B::XYLT1::AP2B1::DMD" into https://www.disgenet.org/search. Output file was originally named "2046__5236__129684__389072__5067__4131__2870__2740__1310__286133__10541__89797__283455__114795__64131__163__1756_gene_gda_summary.tsv" and was downloaded on 9/12/2023. 

 - `disease_gene_associations_PCA_smoking.tsv`: Output from entering "CAMTA1::ECE1::HSPG2::CSMD2::NFIA::LPAR3::CHIA::KCNN3::PAPPA2::SOX13::PTPN14::TP53BP2::DNAH14::PCNXL2::GPR137B::RYR2::OR2L13::GRHL1::ALK::NDUFAF7::GALM::SLC8A1::PRKCE::STON1-GTF2A1L::GPR75-ASB3::EML6::ARHGAP25::LRRTM4::CTNNA2::ACVR1::STK39::CCDC173::PDE11A::CCDC141::PAX3::SRGAP3::DYNC1LI1::ROBO2::GUCA1C::IGSF11::ADCY5::NEK11::MLF1::PHC3::TNIK::KCNMB2::EIF2B5::EIF2B5::LPP::DLG1::CLNK::KCNIP4::GPR125::GABRB1::FIP1L1::EPHA5::DCLK2::KIAA0922::RAPGEF2::TRIO::MYO10::EGFLAM::PDE4D::ARHGEF28::ARSB::EFNA5::GRAMD3::UBE2D2::STK32A::COL23A1::F13A1::DEK::CDKAL1::OPN5::EYS::EPHA7::PREP::CLVS2::TRDN::PERP::HIVEP2::SYNE1::SMOC2::GET4::TNRC18::HDAC9::JAZF1::AMPH::GLI3::WBSCR17::GTPBP10::COL26A1::TBXAS1::CNTNAP2::CSMD1::SGCZ::CSGALNACT1::PIWIL2::DPYSL2::ELP3::ZMAT4::XKR4::CNGB3::SDC2::NCALD::MTSS1::ASAP1::KCNQ3::SLA::ZFAT::KANK1::AK3::OSTF1::PCSK5::GNAQ::TLE1::STX17::ZNF462::SUSD1::COL27A1::AKNA::BRINP1::DENND1A::CAMK1D::C1QL3::RSU1::PARD3::TSPAN14::HPS1::CTBP2::C10orf90::TMEM41B::INSC::NELL1::SLC17A6::SHANK2::DLG2::RAB38::OPCML::RAD52::SLC2A14::GRIN2B::TMCC3::ANKS1B::NUP37::BRAP::DTX1::CCDC92::TMEM132D::CRYL1::ATP8A2::FREM2::FAM124A::PCDH9::ABCC4::NALCN::COL4A1::DCUN1D2::NPAS3::EGLN3::RGS6::CDC42BPB::ATP10A::GABRB3::GABRG3::MEIS2::FBN1::CGNL1::RORA::VPS13C::ZNF609::KIAA1024::PDE8A::MFGE8::ABHD2::CERS3::ADAMTS18::WWOX::CDH13::STX8::PRKCA::SLC39A11::RIT2::ST8SIA5::DCC::MBP::ATP9B::ATCAY::TICAM1::C19orf47::PSG8::PLAUR::ZNF227::VASP::SLC17A7::ZNF175::ZBTB45::PTPRT::CDH22::RUNX1::YBEY::CECR2::CABIN1::SGSM1::LARGE::MPPED1::PJA1" into https://www.disgenet.org/search. Output file was originally named "23261__1889__3339__114784__4774__23566__27159__3782__60676__9580__5784__7159__127602__7107__6262__284521__29841__238__55471__130589__6546__5581__286749__100302652__400954__9938__80059__1496__90__27347__50940__2.tsv" and was downloaded on 9/12/2023. 

 - `miEAA.csv`: Output from entering `MIRS_smoking.txt`into https://ccb-compute2.cs.uni-saarland.de/mieaa2/user_input/ by following manual steps 1 and 2 in `step11b_get_enrichment_all.R`. Was originally named "miEAA - miRNA Enrichment and Annotation -- Analysis results.csv" and was downloaded on 5/29/2023.  

 - `rsIDs_GxAlcohol_effects.txt`: FUMA input. Refer to `Important Note on FUMA details`. 

 - `rsIDs_GxAlcohol_effects_pvals.txt`: FUMA input. Refer to `Important Note on FUMA details`.

 - `rsIDs_GxExercise_effects.txt`: Not used as FUMA input due to there being only one SNP hit. 

 - `rsIDs_GxExercise_effects_pvals.txt`: Not used as FUMA input due to there being only one SNP hit. 

 - `rsIDs_GxGender_effects.txt`: FUMA input. Refer to `Important Note on FUMA details`.

 - `rsIDs_GxGender_effects_pvals.txt`: FUMA input. Refer to `Important Note on FUMA details`.

 - `rsIDs_GxSmoking_effects.txt`: FUMA input. Refer to `Important Note on FUMA details`.

 - `rsIDs_GxSmoking_effects_pvals.txt`: FUMA input. Refer to `Important Note on FUMA details`.

 - `rsIDs_main_effects.txt`: FUMA input. Refer to `Important Note on FUMA details`.

 - `rsIDs_main_effects_pvals.txt`: FUMA input. Refer to `Important Note on FUMA details`.

 - `SNP_MAFs_rsIDs_GxAlcohol_effects.txt`: for possible future use

 - `SNP_MAFs_rsIDs_GxExercise_effects.txt`: for possible future use

 - `SNP_MAFs_rsIDs_GxGender_effects.txt`: for possible future use

 - `SNP_MAFs_rsIDs_GxSmoking_effects.txt`: used by `step11f_miRNA_enrichment_analysis.R` and for possible future use

 - `SNP_MAFs_rsIDs_main_effects.txt`: for possible future use

 - `step11e_logistic_PCA_GxAlcohol_rsIDs_known.txt`: Shows SNP hit subsets used by `step10e_access_common_SNPs.py` to make table 1a.

 - `step11e_logistic_PCA_GxAlcohol_rsIDs_novel.txt`: Shows SNP hit subsets used by `step10e_access_common_SNPs.py` to make table 1a.

 - `step11e_logistic_PCA_GxGender_rsIDs_known.txt`: Shows SNP hit subsets used by `step10e_access_common_SNPs.py` to make table 1a.

 - `step11e_logistic_PCA_GxGender_rsIDs_novel.txt`: Shows SNP hit subsets used by `step10e_access_common_SNPs.py` to make table 1a.

 - `step11e_logistic_PCA_GxSmoking_rsIDs_known.txt`: Shows SNP hit subsets used by `step10e_access_common_SNPs.py` to make table 1a. Used by `step11f_miRNA_enrichment_analysis.R` to make input for https://www.disgenet.org/search.

 - `step11e_logistic_PCA_GxSmoking_rsIDs_novel.txt`: Shows SNP hit subsets used by `step10e_access_common_SNPs.py` to make table 1a. Used by `step11f_miRNA_enrichment_analysis.R` to make input for https://www.disgenet.org/search.

 - `step11e_logistic_PCA_main_rsIDs_known.txt`: Shows SNP hit subsets used by `step10e_access_common_SNPs.py` to make table 1a.

 - `step11e_logistic_PCA_main_rsIDs_novel.txt`: Shows SNP hit subsets used by `step10e_access_common_SNPs.py` to make table 1a.

 - `step11e_NN_GxAlcohol_rsIDs_known.txt`: Shows SNP hit subsets used by `step10e_access_common_SNPs.py` to make table 1a.

 - `step11e_NN_GxAlcohol_rsIDs_novel.txt`: Shows SNP hit subsets used by `step10e_access_common_SNPs.py` to make table 1a.

 - `step11e_NN_GxGender_rsIDs_known.txt`: Shows SNP hit subsets used by `step10e_access_common_SNPs.py` to make table 1a.

 - `step11e_NN_GxGender_rsIDs_novel.txt`: Shows SNP hit subsets used by `step10e_access_common_SNPs.py` to make table 1a.

 - `step11e_NN_GxSmoking_rsIDs_known.txt`: Shows SNP hit subsets used by `step10e_access_common_SNPs.py` to make table 1a. Used by `step11f_miRNA_enrichment_analysis.R` to make input for https://www.disgenet.org/search.

 - `step11e_NN_GxSmoking_rsIDs_novel.txt`: Shows SNP hit subsets used by `step10e_access_common_SNPs.py` to make table 1a. Used by `step11f_miRNA_enrichment_analysis.R` to make input for https://www.disgenet.org/search.

 - `step11e_NN_main_rsIDs_known.txt`: Shows SNP hit subsets used by `step10e_access_common_SNPs.py` to make table 1a.

 - `step11e_NN_main_rsIDs_novel.txt`: Shows SNP hit subsets used by `step10e_access_common_SNPs.py` to make table 1a.

 - `step11e_PCA_GxAlcohol_rsIDs_known.txt`: Shows SNP hit subsets used by `step10e_access_common_SNPs.py` to make table 1a.

 - `step11e_PCA_GxAlcohol_rsIDs_novel.txt`: Shows SNP hit subsets used by `step10e_access_common_SNPs.py` to make table 1a.

 - `step11e_PCA_GxGender_rsIDs_known.txt`: Shows SNP hit subsets used by `step10e_access_common_SNPs.py` to make table 1a.

 - `step11e_PCA_GxGender_rsIDs_novel.txt`: Shows SNP hit subsets used by `step10e_access_common_SNPs.py` to make table 1a.

 - `step11e_PCA_GxSmoking_rsIDs_known.txt`: Shows SNP hit subsets used by `step10e_access_common_SNPs.py` to make table 1a. Used by `step11f_miRNA_enrichment_analysis.R` to make input for https://www.disgenet.org/search.

 - `step11e_PCA_GxSmoking_rsIDs_novel.txt`: Shows SNP hit subsets used by `step10e_access_common_SNPs.py` to make table 1a. Used by `step11f_miRNA_enrichment_analysis.R` to make input for https://www.disgenet.org/search.

 - `step11e_PCA_main_rsIDs_known.txt`: Shows SNP hit subsets used by `step10e_access_common_SNPs.py` to make table 1a.

 - `step11e_PCA_main_rsIDs_novel.txt`: Shows SNP hit subsets used by `step10e_access_common_SNPs.py` to make table 1a.

### Output Files

 - `step11f_logistic_PCA_miRNA_associated_genic_SNP_enrichment.txt`: counts and p-value for enrichment of logistic PCA SNP hits inside of genes linked to miRNA

 - `step11f_NN_miRNA_associated_genic_SNP_enrichment.txt`: counts and p-value for enrichment of NN SNP hits inside of genes linked to miRNA

 - `step11f_PCA_miRNA_associated_genic_SNP_enrichment.txt`: counts and p-value for enrichment of PCA SNP hits inside of genes linked to miRNA

 - `table1a.txt`: refer to manuscript

 - `table1b.txt`: refer to manuscript

 - `table1c.txt`: refer to manuscript

 - `table1d.txt`: refer to manuscript

 - `table2a.txt`: refer to manuscript

 - `table2b.txt`: refer to manuscript

 - `table2c.txt`: refer to manuscript

 - `table3a.txt`: refer to manuscript

 - `table3b.txt`: refer to manuscript

 - `table3c.txt`: refer to manuscript

 - `table3d.txt`: refer to manuscript

### Resource Folders

 - `all_sig_rsIDs_logistic_PCA`: all statistically significant SNP hits for logistic PCA latent phenotypes

 - `all_sig_rsIDs_NN`: all statistically significant SNP hits for NN latent phenotypes

 - `all_sig_rsIDs_PCA`: all statistically significant SNP hits for PCA latent phenotypes

 - `genevestigator_hits`: all studies with differentially expressed genes

 - `genevestigator_hits_methods`: contains methodological details regarding `genevestigator_hits`

