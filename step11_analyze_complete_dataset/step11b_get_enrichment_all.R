#Turns out you need to install bioconductor related things this way.
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.16")

#BiocManager::install("limma")
#BiocManager::install("AnnotationDbi")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("msigdbr")
#BiocManager::install("GO.db")
#BiocManager::install("clusterProfiler")
#install.packages("tidyr")
#install.packages("gt")
#install.packages("webshot2")
#install.packages("tidytext")
#install.packages("dplyr")
#install.packages("VennDiagram")
#install.packages("stringr") 

packs = c("ggplot2", "limma", "tidyr", "AnnotationDbi", "org.Hs.eg.db", 
          "GO.db", "msigdbr", "clusterProfiler", "gt", "webshot2", "tidytext", 
          "dplyr", "VennDiagram", "stringr")
invisible(lapply(packs, library, character.only = TRUE))

names1 = c("GxSmoking", "GxAlcohol", "GxGender", "main")
names2 = c("smoking", "alcohol", "gender", "main")
hit_gene_pair_sets = list()
gene_ID_dicts = list()
info_KEGG = list()
info_GO = list()
info_REG = list()
for (i in 1:4){
    
    fname1 = paste("rsIDs_", names1[i], "_effects.txt", sep = "")
    fname2 = paste("annov_", names2[i], ".txt", sep = "")
    fname3 = paste("gene_list_", names2[i], ".txt", sep = "")
    my_SNPs <- read.csv(fname1, header = T, sep = "\t")
    my_SNP_gene_pairs <- read.csv(fname2, header = T, sep = "\t") %>% 
      tidyr::separate("uniqID", c("chr", "pos", NA, NA)) %>% 
      merge(my_SNPs, by = c("chr", "pos"))
    hit_gene_pair_sets[[i]] <- my_SNP_gene_pairs
    write.table(my_SNP_gene_pairs[, "symbol"], fname3, row.names = F, col.names = F)

    my_gene_IDs <- my_SNP_gene_pairs[my_SNP_gene_pairs["dist"] < 300000, ]$gene
    Enterez_IDs <- mapIds(org.Hs.eg.db, keys = my_gene_IDs, keytype="ENSEMBL", column = "ENTREZID")
    my_Enterez_IDs <- Enterez_IDs[is.na(Enterez_IDs) == F]
    gene_ID_dicts[[i]] <- my_Enterez_IDs

    KEGG_genes <- msigdbr(species = "human", subcategory = "KEGG")
    KEGG_t2g <- KEGG_genes %>% dplyr::distinct(gs_name, entrez_gene) %>% as.data.frame()
    out_KEGG <- enricher(gene = unique(my_Enterez_IDs), TERM2GENE = KEGG_t2g)
    df_KEGG <- out_KEGG@result[out_KEGG@result[6] < 0.05, c(1, 6, 9, 8)]
    if (nrow(df_KEGG) > 0){
        effect_name = paste("Enriched with gene by", names2[i], "effects", sep = " ")
        if (names2[i] == "main"){effect_name = "Enriched with main effects"}
        df_KEGG$effect <- rep(effect_name, nrow(df_KEGG))  
        rownames(df_KEGG) <- NULL
        info_KEGG[[i]] <- df_KEGG}

    GO_genes <- msigdbr(species = "human", category = "C5")
    GO_t2g <- GO_genes %>% dplyr::distinct(gs_name, entrez_gene) %>% as.data.frame()
    out_GO <- enricher(gene = unique(my_Enterez_IDs), TERM2GENE = GO_t2g)
    df_GO <- out_GO@result[out_GO@result[6] < 0.05, c(1, 6, 9, 8)]
    if (nrow(df_GO) > 0){
        effect_name = paste("Enriched with gene by", names2[i], "effects", sep = " ")
        if (names2[i] == "main"){effect_name = "Enriched with main effects"}
        df_GO["effect"] <- rep(effect_name, nrow(df_GO))
        rownames(df_GO) <- NULL
        info_GO[[i]] <- df_GO}

    REG_genes <- msigdbr(species = "human", category = "C3")
    REG_t2g <- REG_genes %>% dplyr::distinct(gs_name, entrez_gene) %>% as.data.frame()
    out_REG <- enricher(gene = unique(my_Enterez_IDs), TERM2GENE = REG_t2g)
    df_REG <- out_REG@result[out_REG@result[6] < 0.05, c(1, 6, 9, 8)]
    if (nrow(df_REG) > 0){
        effect_name = paste("Enriched with gene by", names2[i], "effects", sep = " ")
        if (names2[i] == "main"){effect_name = "Enriched with main effects"}
        df_REG["effect"] <- rep(effect_name, nrow(df_REG))
        rownames(df_REG) <- NULL
        info_REG[[i]] <- df_REG}
}

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
  
  # note at this point I paste the gene symbols together and input the gene
  # input the gene list into the database (https://www.disgenet.org/dbinfo)
  # paste(rsID_data[["symbol"]], collapse = "::")
  # then download the summary file
  
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
  novel_rsiD_tables_both[[i]] <- rsID_data[has_both, ]
  novel_rsiD_tables_pheno[[i]] <- rsID_data[has_gda_only, ]
  novel_rsiD_tables_env[[i]] <- rsID_data[has_DE_only, ]
  novel_rsiD_tables_neither[[i]] <- rsID_data[has_neither, ]
  i = i + 1
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


list_PCA <- novel_rsiD_tables[[1]][["symbol"]]
list_logistic_PCA <- novel_rsiD_tables[[2]][["symbol"]]
list_NN <- novel_rsiD_tables[[3]][["symbol"]]
venn.plot <- draw.triple.venn(
  area1 = length(list_PCA),
  area2 = length(list_logistic_PCA),
  area3 = length(list_NN),
  n12 = length(intersect(list_PCA, list_logistic_PCA)),
  n23 = length(intersect(list_logistic_PCA, list_NN)),
  n13 = length(intersect(list_PCA, list_NN)),
  n123 = length(intersect(list_PCA, intersect(list_logistic_PCA, list_NN))),
  category = c("List 1", "List 2", "List 3"),
  fill = c("blue", "green", "yellow")
)
venn.plot

MIRs <- info_REG[[1]][grepl("MIR",  unlist(info_REG[[1]][1])), 1]
MIRs <- gsub("[ATGC]{2,}_", "", MIRs)
MIRs <- gsub("_MIR", "__MIR", MIRs)
MIRs <- unname(unlist(strsplit(MIRs, split = "__")))
MIRs <- tolower(MIRs)
MIRs <- gsub("mir", "hsa-miR-", MIRs)
MIRs <- gsub("_", "-", MIRs)
MIRs <- data.frame(MIRs)
write.table(MIRs, "MIRS_smoking.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote = FALSE)

unknown_MIRs = read.csv("MIRS_smoking_unknown.txt", header = F, sep = "\t")
MIRs2 = data.frame(setdiff(unname(unlist(MIRs)), unname(unlist(unknown_MIRs))))
write.table(MIRs2, "MIRS2_smoking.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote = FALSE)

TFs = info_REG[[3]][grepl("MIR",  unlist(info_REG[[3]][1])) == FALSE, 1]
TFs = TFs[grepl("LET",  TFs) == FALSE]
TFs = sapply(strsplit(TFs, split = "_"), `[`, 1)
TFs_Enterez_IDs <- mapIds(org.Hs.eg.db, keys = TFs, keytype="SYMBOL", column = "ENTREZID")
TFs_Enterez_IDs <- TFs_Enterez_IDs[is.na(TFs_Enterez_IDs) == F]

TFs_KEGG_genes <- msigdbr(species = "human", subcategory = "BP")
TFs_KEGG_t2g <- TFs_KEGG_genes %>% dplyr::distinct(gs_name, entrez_gene) %>% as.data.frame()
TFs_out_KEGG <- enricher(gene = unique(TFs_Enterez_IDs), TERM2GENE = TFs_KEGG_t2g)
TFs_info_KEGG_gender <- TFs_out_KEGG@result[TFs_out_KEGG@result[6] < 0.05, c(1, 6, 9)]

TFs_info_KEGG_gender[,"ID"] <- gsub("GOBP_", "", TFs_info_KEGG_gender[,"ID"])
TFs_info_KEGG_gender[,"ID"] <- gsub("_", " ", TFs_info_KEGG_gender[,"ID"])
TFs_info_KEGG_gender[,"ID"] <- tolower(TFs_info_KEGG_gender[,"ID"])

gender_table <- TFs_info_KEGG_gender[1:10, ]
gender_table <- gt(gender_table, 
                   rowname_col = "ID") %>% 
  cols_label(
    p.adjust = md("Adjusted p-value"),
    Count = md("Number of genes"),
  ) %>%
  fmt_scientific(
    columns = c("p.adjust"),
    decimals = 2
  ) %>% 
  cols_align(
    columns = c("p.adjust"),
    align = "right"
  ) %>% 
  cols_align(
    columns = c("Count"),
    align = "right"
  ) %>% 
  tab_stubhead(label = "GO gene set ID") %>% 
  cols_width(
    ID ~ px(350),
    p.adjust ~ px(175),
    Count ~ px(175),
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
gender_table
gtsave(gender_table, "gender_table.docx")
gtsave(gender_table, "gender_table.png")

alcohol_table <- info_KEGG[[2]][ ,c("ID", "p.adjust", "geneID")]
alcohol_table["geneID"] <- c("ALDH5", "ALDH5", "PIK3C2G", "PIK3C2G", "ID4")
alcohol_table[,"ID"] <- gsub("KEGG_", "", alcohol_table[,"ID"])
alcohol_table[,"ID"] <- gsub("_", " ", alcohol_table[,"ID"])
alcohol_table[,"ID"] <- tolower(alcohol_table[,"ID"])
alcohol_table <- gt(alcohol_table, 
                    rowname_col = "ID") %>% 
  cols_label(
    p.adjust = md("Adjusted p-value"),
    geneID = md("gene name"),
  ) %>%
  fmt_scientific(
    columns = c("p.adjust"),
    decimals = 2
  ) %>% 
  cols_align(
    columns = c("p.adjust"),
    align = "right"
  ) %>% 
  cols_align(
    columns = c("geneID"),
    align = "right"
  ) %>% 
  tab_stubhead(label = "KEGG gene set ID") %>% 
  cols_width(
    ID ~ px(350),
    p.adjust ~ px(175),
    geneID ~ px(175),
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
alcohol_table
gtsave(alcohol_table, "alcohol_table.docx")
gtsave(alcohol_table, "alcohol_table.png")

main_table <- info_GO[[4]][1:10 ,c("ID", "p.adjust", "Count")]
main_table[,"ID"] <- gsub("GO_", "", main_table[,"ID"])
main_table[,"ID"] <- gsub("_", " ", main_table[,"ID"])
main_table[,"ID"] <- tolower(main_table[,"ID"])
main_table <- gt(main_table, 
                 rowname_col = "ID") %>% 
  cols_label(
    p.adjust = md("Adjusted p-value"),
    Count = md("Number of genes"),
  ) %>%
  fmt_scientific(
    columns = c("p.adjust"),
    decimals = 2
  ) %>% 
  cols_align(
    columns = c("p.adjust"),
    align = "right"
  ) %>% 
  cols_align(
    columns = c("Count"),
    align = "right"
  ) %>% 
  tab_stubhead(label = "GO gene set ID") %>% 
  cols_width(
    ID ~ px(350),
    p.adjust ~ px(175),
    Count ~ px(175),
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
main_table
gtsave(main_table, "main_table.docx")
gtsave(main_table, "main_table.png")

miEAA_hits <- read.csv("miEAA.csv", header = T, sep = ",")
smoking_table <- miEAA_hits[1:10, c("Subcategory", "P.adjusted", "Observed")]
smoking_table <- gt(smoking_table, 
                    rowname_col = "Subcategory") %>% 
  cols_label(
    P.adjusted = md("Adjusted p-value"),
    Observed = md("Number of miRNA"),
  ) %>%
  fmt_scientific(
    columns = c("P.adjusted"),
    decimals = 2
  ) %>% 
  cols_align(
    columns = c("P.adjusted"),
    align = "right"
  ) %>% 
  cols_align(
    columns = c("Observed"),
    align = "right"
  ) %>% 
  tab_stubhead(label = "KEGG gene set ID") %>% 
  cols_width(
    Subcategory ~ px(350),
    P.adjusted ~ px(175),
    Observed ~ px(175),
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
smoking_table
gtsave(smoking_table, "smoking_table.docx")
gtsave(smoking_table, "smoking_table.png")

# this was right before the fisher exact test
# This was for checking enrichment of phenotypes for SNPs with heterotic effects
# There appear to be none. 
# list0 <- info_REG[[1]][grepl("MIR",  unlist(info_REG[[1]][1])), c("ID", "geneID")]
# list0_mirs <- list0[grepl("MIR",  unlist(list0["ID"])), "ID"]
# list0_mirs <- gsub("[ATGC]{2,}_", "", list0_mirs)
# list0_mirs <- gsub("_MIR", "__MIR", list0_mirs)
# list0_mirs <- strsplit(list0_mirs, split = "__")
# lengths <- sapply(list0_mirs, length)
# list0_geneIDs_values <- list0[grepl("MIR",  unlist(list0["ID"])), "geneID"]
# list0_geneIDs <- list()
# for (i in 1:length(lengths)){
#   list0_geneIDs[[i]] <- rep(list0_geneIDs_values[i], lengths[i]) 
# }
# list0_mirs <- unlist(list0_mirs)
# list0_mirs <- tolower(list0_mirs)
# list0_mirs <- gsub("mir", "hsa-miR-", list0_mirs)
# list0_mirs <- gsub("_", "-", list0_mirs)
# list0_keepers <- list0_mirs %in% HCM_mirs
# list0 <- data.frame(list0_mirs , unlist(list0_geneIDs))
# colnames(list0) <- c("ID", "geneID")
# list0 <- list0[list0_keepers, "geneID"]

