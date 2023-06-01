# WARNING: this code contains some manual steps

#---------------------------------------------------------------------------------------------------------
#first run this code.
#---------------------------------------------------------------------------------------------------------
MIRs <- info_REG[[1]][grepl("MIR",  unlist(info_REG[[1]][1])), 1]
MIRs <- gsub("[ATGC]{2,}_", "", MIRs)
MIRs <- gsub("_MIR", "__MIR", MIRs)
MIRs <- unname(unlist(strsplit(MIRs, split = "__")))
MIRs <- tolower(MIRs)
MIRs <- gsub("mir", "hsa-miR-", MIRs)
MIRs <- gsub("_", "-", MIRs)
MIRs <- data.frame(MIRs)
write.table(MIRs, "MIRS_smoking.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote = FALSE)
#---------------------------------------------------------------------------------------------------------



#---------------------------------------------------------------------------------------------------------
# MANUAL STEP
# paste the miRNA names from MIRS_smoking.txt into the appropriate step at
# https://ccb-compute2.cs.uni-saarland.de/mieaa2/user_input/
# Then you'll get a list of miRNA names not on record.
# Put that file in your working directory and rename it "MIRS_smoking_unknown.txt"
# Then run this code
#---------------------------------------------------------------------------------------------------------
unknown_MIRs = read.csv("MIRS_smoking_unknown.txt", header = F, sep = "\t")
MIRs2 = data.frame(setdiff(unname(unlist(MIRs)), unname(unlist(unknown_MIRs))))
write.table(MIRs2, "MIRS2_smoking.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote = FALSE)

HCM_mirs <- read.csv("miEAA.csv", header = T, sep = ",")
name = HCM_mirs[9, "Subcategory"] 
row <- (HCM_mirs["Subcategory"] == name)
HCM_mirs <- unlist(strsplit(HCM_mirs[row, "miRNAs.precursors"], "; "))
list0 <-info_REG[[1]][grepl("MIR",  unlist(info_REG[[1]][1])), "geneID"]

#input: the list of all gene numbers associated with at least 1 miRNA
#output the list of all unique rsIDs/gen loc pairs near genes that enrich miRNA gene sets
mir_gene_numbers <- names(table(unlist(strsplit(list0, "/"))))
ID_dict <- gene_ID_dicts[[1]]
mir_gene_ensgIDs <- names(ID_dict[match(mir_gene_numbers, ID_dict)])
all_gene_ensgIDs <- unlist(hit_gene_pair_sets[[1]]["gene"])
mir_inds_in_all <- match(mir_gene_ensgIDs, all_gene_ensgIDs)
cols <- c("rsID", "chr", "annot")
all_rsids <- hit_gene_pair_sets[[1]][cols]
mir_rsids <-  all_rsids[mir_inds_in_all, cols]
all_rsids <- all_rsids[duplicated(all_rsids) == FALSE, ]
mir_rsids <- mir_rsids[duplicated(mir_rsids) == FALSE, ]
write.table(mir_rsids, "MIRS_rsids.txt", row.names = F, col.names = T, sep="\t", quote = F)
#---------------------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------------------
# MANUAL STEP
# run this code to print 3 lists of :: seperated gene symbols.
# the 1st, 2nd, and 3rd lists are for PCA, logistic PCA, and NN respectively.
# input each gene list into the database (https://www.disgenet.org/dbinfo)
# make sure "genes" is selected in the options above the search box
# then download each summary of gene-disease associations file.
# name them  "disease_gene_associations_PCA_smoking.tsv",
#            "disease_gene_associations_logistic_PCA_smoking.tsv", and
#            "disease_gene_associations_NN_smoking.tsv" respectively. 
#---------------------------------------------------------------------------------------------------------
for (model in c("PCA", "logistic_PCA", "NN")) {
  fname <- paste("step11e_", model, "_GxSmoking_rsIDs_known.txt", sep = "")
  rsIDs_model_known <- read.csv(fname, header = T, sep = "\t")
  fname <- paste("step11e_", model, "_GxSmoking_rsIDs_novel.txt", sep = "")
  rsIDs_model_novel <- read.csv(fname, header = T, sep = "\t")
  rsIDs_model <- rbind(rsIDs_model_known, rsIDs_model_novel)
  keepers <- match(rsIDs_model[["rsID"]], mir_rsids[["rsID"]])
  keepers <- keepers[is.na(keepers) == FALSE]
  mir_rsids_i <- mir_rsids[keepers, ]
  
  fig_rsIDs <- mir_rsids_i[mir_rsids_i[, "annot"] == "intronic", "rsID"]
  fig_keepers <- match(fig_rsIDs, rsIDs_model_novel[["rsID"]])
  fig_keepers <- fig_keepers[is.na(fig_keepers) == FALSE]
  fig_rsIDs <- rsIDs_model_novel[fig_keepers, ]
  rsID_data <- hit_gene_pair_sets[[1]]
  data_keepers <- match(fig_rsIDs, rsID_data[["rsID"]])
  data_keepers <- data_keepers[is.na(data_keepers) == FALSE]
  rsID_data <- rsID_data[data_keepers, c("rsID", "chr", "pos", "symbol", "annot")]
  rsID_data <- rsID_data[rsID_data[["annot"]] == "intronic", c("rsID", "chr", "pos", "symbol")]
  print(paste(rsID_data[["symbol"]], collapse = "::"))
}
#---------------------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------------------
# MANUAL STEP
# First use genevestigator to get the files in genevestigator_hits
# then run the code below
#---------------------------------------------------------------------------------------------------------
DE_fnames_all <- list.files(path = "genevestigator_hits", full.names = T)
DE_fnames_status <- DE_fnames_all[grepl("smoker_status", DE_fnames_all)]
DE_fnames_exposure <- DE_fnames_all[grepl("smoke_exposure", DE_fnames_all)]
DE_fname_sets <- list(DE_fnames_all, DE_fnames_status, DE_fnames_exposure)
setnames <- c("all", "smoker_status", "smoke_exposure")
DE_data_sets = list()
for (i in 1:3){
  DE_fnames <- DE_fname_sets[[i]]
  DE_files <- lapply(DE_fnames, function(name){
    file_data <- read.csv(name, header = F, sep = "\t")
    file_data$name <- name
    return(file_data)})
  raw_DE_data <- Reduce(rbind, DE_files)
  DE_data <- Reduce(rbind, DE_files) %>%
    group_by(V1) %>%
    summarise(average = mean(V2),
              average_abs = mean(abs(V2)),
              stdev = sd(V2),
              num_hits = n()) %>%
    rename(symbol = V1)
  if (setnames[i] == "all"){
    keepers1 <- DE_data["num_hits"] >= 2 & DE_data["average_abs"] >= 0.6
    keepers2 <- DE_data["num_hits"] >= 3
    keepers <- keepers1 | keepers2
    keeper_genes <- DE_data$symbol[keepers]}
  else{
    keepers <- DE_data$symbol %in% keeper_genes
  }
  colnames(DE_data)[2] <- paste(colnames(DE_data)[2], setnames[i], sep = "_")
  colnames(DE_data)[3] <- paste(colnames(DE_data)[3], setnames[i], sep = "_")
  colnames(DE_data)[4] <- paste(colnames(DE_data)[4], setnames[i], sep = "_")
  colnames(DE_data)[5] <- paste(colnames(DE_data)[5], setnames[i], sep = "_")
  DE_data_sets[[i]] <- DE_data[keepers, ]
}
#---------------------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------------------
# run the code below to get the fishers tables for miRNA SNP location analysis
# also to get the analysis of known importance
#---------------------------------------------------------------------------------------------------------
tests <- list()
tables <- list()
novel_rsiD_tables_both <- list()
novel_rsiD_tables_pheno <- list()
novel_rsiD_tables_env <- list()
novel_rsiD_tables_neither <- list()
i = 1
for (model in c("PCA", "logistic_PCA", "NN")) {
  fname <- paste("step11e_", model, "_GxSmoking_rsIDs_known.txt", sep = "")
  rsIDs_model_known <- read.csv(fname, header = T, sep = "\t")
  fname <- paste("step11e_", model, "_GxSmoking_rsIDs_novel.txt", sep = "")
  rsIDs_model_novel <- read.csv(fname, header = T, sep = "\t")
  rsIDs_model <- rbind(rsIDs_model_known, rsIDs_model_novel)
  keepers <- match(rsIDs_model[["rsID"]], mir_rsids[["rsID"]])
  keepers <- keepers[is.na(keepers) == FALSE]
  mir_rsids_i <- mir_rsids[keepers, ]
  keepers <- match(rsIDs_model[["rsID"]], all_rsids[["rsID"]])
  keepers <- keepers[is.na(keepers) == FALSE]
  all_rsids_i <- all_rsids[keepers, ]
  
  mir_genic1 <- sum(mir_rsids_i[, "annot"] == "intronic")
  mir_genic2 <- sum(mir_rsids_i[, "annot"] == "exonic")
  mir_genic3 <- sum(mir_rsids_i[, "annot"] == "UTR5")
  mir_genic4 <- sum(mir_rsids_i[, "annot"] == "UTR3")
  mir_genic <- mir_genic1 + mir_genic2 + mir_genic3 + mir_genic4
  mir_not_genic <- length(mir_rsids_i[, "annot"]) - mir_genic
  all_genic1 <- sum(all_rsids_i[, "annot"] == "intronic")
  all_genic2 <- sum(all_rsids_i[, "annot"] == "exonic")
  all_genic3 <- sum(all_rsids_i[, "annot"] == "UTR5")
  all_genic4 <- sum(all_rsids_i[, "annot"] == "UTR3")
  all_genic <- all_genic1 + all_genic2 + all_genic3 + all_genic4
  all_not_genic <- length(all_rsids_i[, "annot"]) - all_genic
  not_mir_genic <- all_genic - mir_genic
  not_mir_not_genic <- all_not_genic - mir_not_genic
  mir <- c(mir_genic, mir_not_genic)
  not_mir <- c(not_mir_genic, not_mir_not_genic)
  data <- data.frame(mir, not_mir, row.names = c("genic", "not_genic"))
  colnames(data) <- c("mir", "not_mir")
  tables[[i]] <- data
  tests[[i]] <- fisher.test(data)
  
  fig_rsIDs <- mir_rsids_i[mir_rsids_i[, "annot"] == "intronic", "rsID"]
  fig_keepers <- match(fig_rsIDs, rsIDs_model_novel[["rsID"]])
  fig_keepers <- fig_keepers[is.na(fig_keepers) == FALSE]
  fig_rsIDs <- rsIDs_model_novel[fig_keepers, ]
  rsID_data <- hit_gene_pair_sets[[1]]
  data_keepers <- match(fig_rsIDs, rsID_data[["rsID"]])
  data_keepers <- data_keepers[is.na(data_keepers) == FALSE]
  rsID_data <- rsID_data[data_keepers, c("rsID", "chr", "pos", "symbol", "annot")]
  rsID_data <- rsID_data[rsID_data[["annot"]] == "intronic", c("rsID", "chr", "pos", "symbol")]
  rownames(rsID_data) <- 1:nrow(rsID_data)
  
  fname <- paste("disease_gene_associations_", model, "_smoking.tsv", sep = "")
  diseases <- read.csv(fname, header = T, sep = "\t")
  CV_inds <- grep("Cardiovascular Diseases", diseases[["Disease_Class"]])
  diseases <- diseases[CV_inds, ]
  disease_genes <- diseases[, c("Gene", "Score_gda", "Disease")]
  disease_genes <- disease_genes %>% 
    group_by(Gene) %>% 
    summarize(gda_score_sum = sum(Score_gda),
              diseases = paste(tolower(Disease), collapse = "::")) %>%
    mutate(diseases = sapply(str_split(diseases, "::"), 
                             function(x) paste(x[1:min(10, length(x))], 
                                               collapse = "; ")
    )
    ) %>%
    as.data.frame
  
  fname <- paste("all_sig_rsIDs_", model, "/rsIDs_GxSmoking_effects.txt", sep = "")
  original_snps <- read.csv(fname, header = T, sep = "\t")
  fname <- "SNP_MAFs_rsIDs_GxSmoking_effects.txt"
  snp_mafs <- read.csv(fname, header = T, sep = "\t")[c("rsID", "MAF")]
  snp_mafs <- snp_mafs[duplicated(snp_mafs) == FALSE,]
  snp_mafs["MAF"] <- sapply(snp_mafs[["MAF"]], function(x) min(x, 1 - x))
  snp_info <- merge(original_snps, snp_mafs, by = "rsID")
  snp_info <- snp_info %>% 
    group_by(rsID) %>% 
    summarize(best_p_value = max(pEDGE2), 
              maf = max(MAF)) %>% #note: SNPs with the same rsID have the same maf.
    as.data.frame
  
  colnames(disease_genes) <- c("symbol", "gda_score_sum", "diseases")
  rsID_data <- merge(rsID_data, snp_info, by = "rsID")
  rsID_data <- merge(rsID_data, disease_genes, by = "symbol", all.x = TRUE)
  rsID_data <- rsID_data[order(desc(rsID_data[, "gda_score_sum"])), ]
  for (j in 1:3){
    DE_data <- DE_data_sets[[j]]
    rsID_data <- merge(rsID_data, DE_data, by = "symbol", all.x = TRUE)
  }
  sorted_indices <- order(rsID_data$gda_score_sum,
                          rsID_data$num_hits_all, 
                          rsID_data$average_abs_all,
                          rsID_data$best_p_value,
                          decreasing = TRUE)
  rsID_data <- rsID_data[sorted_indices, ]
  rownames(rsID_data) <- as.character(1:nrow(rsID_data))
  
  has_gda_sum <- is.na(rsID_data$gda_score_sum) == FALSE
  has_DE_effect <- is.na(rsID_data$num_hits_all) == FALSE
  has_both <- has_gda_sum & has_DE_effect
  has_gda_only <- has_gda_sum & (has_DE_effect == FALSE)
  has_DE_only <- (has_gda_sum == FALSE) & has_DE_effect
  has_neither <- (has_gda_sum == FALSE) & (has_DE_effect == FALSE)
  rsID_data["phenotype model"] <- model
  novel_rsiD_tables_both[[i]] <- rsID_data[has_both, ]
  novel_rsiD_tables_pheno[[i]] <- rsID_data[has_gda_only, ]
  novel_rsiD_tables_env[[i]] <- rsID_data[has_DE_only, ]
  novel_rsiD_tables_neither[[i]] <- rsID_data[has_neither, ]
  i = i + 1
}

both_top_five <- Reduce(rbind, novel_rsiD_tables_both) %>% 
                 filter(.data$average_abs_all > 0.4) %>% 
                 arrange(desc(gda_score_sum)) %>% 
                 select("symbol", "rsID", "phenotype model",  
                        "best_p_value", "maf", "gda_score_sum", 
                        "num_hits_all", "average_all") %>% 
                 rename(`gene symbol` = symbol, `best p value` = best_p_value,
                       `gda score sum` = gda_score_sum,
                       `average DE log2(fold change)` = average_all, 
                       `number of DE reps` = num_hits_all) %>% 
                 slice(1:5)
both_roles <- c("apoptosis",
                "G protein",
                "CVD-related mRNA regulation",
                "CVD-related transcription regulation",
                "multiple CVD-related pathways")

pheno_top_five <- Reduce(rbind, novel_rsiD_tables_pheno) %>% 
                  arrange(desc(gda_score_sum)) %>% 
                  select("symbol", "rsID", "phenotype model", 
                         "best_p_value", "maf", "gda_score_sum") %>% 
                  rename(`gene symbol` = symbol, 
                         `best p value` = best_p_value,
                         `gda score sum` = gda_score_sum) %>% 
                  slice(1:5)
pheno_roles <- c("cardiac calcium channel",
                 "estrogen receptor",
                 "endothelin converting enzyme",
                 "cardiac muscle contraction",
                 "perlecan protein")

env_top_five <- Reduce(rbind, novel_rsiD_tables_env) %>% 
                filter(.data$average_abs_all > 0.4) %>% 
                arrange(desc(num_hits_all), 
                        desc(average_abs_all)) %>% 
                select("symbol", "rsID", "phenotype model", 
                       "best_p_value", "maf", 
                       "num_hits_all", "average_all") %>% 
                rename(`gene symbol` = symbol, 
                       `best p value` = best_p_value,
                       `average DE log2(fold change)` = average_all, 
                       `number of DE reps` = num_hits_all) %>% 
                slice(1:5)
env_roles  <- c("neurogenesis",
                "chondroitin sulfate assembly",
                "cytoskeleton formation",
                "chloride transport",
                "potassium channels")


neither_top_5 <- Reduce(rbind, novel_rsiD_tables_neither) %>% 
                 arrange(best_p_value) %>%  
                 select("symbol", "rsID", "phenotype model", 
                        "best_p_value", "maf") %>% 
                 rename(`gene symbol` = symbol, 
                        `best p value` = best_p_value) %>% 
                 slice(1:5)
neither_roles <- c("Vesicle trafficking protein",
                   "Transcription factor",
                   "glycosylation",
                   "Ras regulator",
                   "mTOR regulator")

gene_functions <- list(both_roles, pheno_roles, env_roles, neither_roles)
tables <- list(both_top_five, pheno_top_five, env_top_five, neither_top_5)
out_tables <- list()
docx_names <- c("top_five_both.docx", "top_five_pheno.docx", 
                "top_five_env.docx", "top_5_neither.docx")
png_names <- c("top_five_both.png", "top_five_pheno.png", 
               "top_five_env.png", "top_5_neither.png")

for (i in 1:4) {

      role <- gene_functions[[i]]
      table <- tables[[i]] %>% 
      add_column(`function` = role, .after = 1)
      docx_name <- docx_names[i]
      png_names <- png_names[i]
      table <- gt(table, 
                  rowname_col = "gene symbol") %>%
      fmt_scientific(
        columns = "best p value",
        decimals = 2
      ) %>% 
      fmt_number(
        columns = "maf",
        decimals = 2
      ) %>% 
      cols_align(
        columns = c("rsID"),
        align = "right"
      ) %>% 
      cols_align(
        columns = c("phenotype model"),
        align = "right"
      ) %>% 
      cols_align(
        columns = c("best p value"),
        align = "right"
      ) %>% 
      cols_align(
        columns = c("maf"),
        align = "right"
      ) %>% 
      tab_stubhead(label = "gene symbol") %>% 
      cols_width(
        `gene symbol` ~ px(125),
        `function` ~ px(200),
        `rsID` ~ px(130),
        `phenotype model` ~ px(125),
        `best p value` ~ px(125),
        `maf` ~ px(75),
      ) %>%
      tab_options(
        column_labels.background.color = "#0b1d78",
        table.background.color = "#d8f9ff") %>%
      tab_style(
        style = cell_text(align = "right", v_align = "middle"),
        locations = cells_stub()
      ) %>%
      tab_style(
        style = cell_text(align = "right", v_align = "middle"),
        locations = cells_stubhead()
      ) %>%
      tab_style(
        style = cell_text(align = "left", v_align = "middle"),
        locations = cells_column_labels()
      ) %>%
      tab_style(
        style = cell_text(align = "left", v_align = "middle"),
        locations = cells_body()
      ) %>%
      opt_horizontal_padding(scale = 3)
    
    if (i %in% c(1,3)){
      table <- fmt_number(table,
      columns = "average DE log2(fold change)",
      decimals = 2
      ) %>%  
      cols_align(
        columns = "average DE log2(fold change)",
        align = "right"
      ) %>% 
      cols_align(
        columns = "number of DE reps",
        align = "right"
      ) %>% 
      cols_width(
        `average DE log2(fold change)` ~ px(160),
        `number of DE reps` ~ px(110)
      )
    }
      
    if (i %in% c(1,2)){
      table <- cols_align(table,
        columns = c("gda score sum"),
        align = "right"
      ) %>% 
      cols_width(
        `gda score sum` ~ px(100),
      )
    }
      
    out_tables[[i]] <- table
    #gtsave(table, docx_name)
    #gtsave(table, png_name)
}

write.table(novel_rsiD_tables_both[[1]], "examples_PCA_smoking_both.txt", row.names = F, col.names = T, sep="\t", quote = F)
write.table(novel_rsiD_tables_pheno[[1]], "examples_PCA_smoking_pheno.txt", row.names = F, col.names = T, sep="\t", quote = F)
write.table(novel_rsiD_tables_env[[1]], "examples_PCA_smoking_env.txt", row.names = F, col.names = T, sep="\t", quote = F)
write.table(novel_rsiD_tables_neither[[1]], "examples_PCA_smoking_neither.txt", row.names = F, col.names = T, sep="\t", quote = F)
write.table(novel_rsiD_tables_both[[2]], "examples_logistic_PCA_smoking_both.txt", row.names = F, col.names = T, sep="\t", quote = F)
write.table(novel_rsiD_tables_pheno[[2]], "examples_logistic_PCA_smoking_pheno.txt", row.names = F, col.names = T, sep="\t", quote = F)
write.table(novel_rsiD_tables_env[[2]], "examples_logistic_PCA_smoking_env.txt", row.names = F, col.names = T, sep="\t", quote = F)
write.table(novel_rsiD_tables_neither[[2]], "examples_logistic_PCA_smoking_neither.txt", row.names = F, col.names = T, sep="\t", quote = F)
write.table(novel_rsiD_tables_both[[3]], "examples_NN_smoking_both.txt", row.names = F, col.names = T, sep="\t", quote = F)
write.table(novel_rsiD_tables_pheno[[3]], "examples_NN_smoking_pheno.txt", row.names = F, col.names = T, sep="\t", quote = F)
write.table(novel_rsiD_tables_env[[3]], "examples_NN_smoking_env.txt", row.names = F, col.names = T, sep="\t", quote = F)
write.table(novel_rsiD_tables_neither[[3]], "examples_NN_smoking_neither.txt", row.names = F, col.names = T, sep="\t", quote = F)
#---------------------------------------------------------------------------------------------------------
