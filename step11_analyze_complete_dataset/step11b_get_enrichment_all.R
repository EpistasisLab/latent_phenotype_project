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

