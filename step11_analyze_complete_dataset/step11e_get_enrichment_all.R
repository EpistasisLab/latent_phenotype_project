#---------------------------------------------------------------------------------------------------------
# WARNING: step11b_get_enrichment_all.R contains some manual steps
# DO NOT run all of the code at once
# Run it all up to each manual step. Then do the manual step.
# Repeat until there are no remaining intermediate manual steps.
#---------------------------------------------------------------------------------------------------------

#Turns out you need to install bioconductor related things this way.
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.16")
# 
# BiocManager::install("limma")
# BiocManager::install("AnnotationDbi")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("msigdbr")
# BiocManager::install("GO.db")
# BiocManager::install("clusterProfiler")
# install.packages("tidyr")
# install.packages("gt")
# install.packages("webshot2")
# install.packages("tidytext")
# install.packages("dplyr")
# install.packages("VennDiagram")
# install.packages("stringr")
# install.packages("MASS")
# install.packages("ggforce")
# install.packages("DiagrammeR")
# install.packages("rsvg")
# install.packages('DiagrammeRsvg')

#---------------------------------------------------------------------------------------------------------
# WARNING: step11b_get_enrichment_all.R contains some manual steps
# DO NOT run all of the code at once
# Run it all up to each manual step. Then do the manual step.
# Repeat until there are no remaining intermediate manual steps.
#---------------------------------------------------------------------------------------------------------

packs = c("ggplot2", "limma", "tidyr", "AnnotationDbi", "org.Hs.eg.db", "DiagrammeR", 
          "GO.db", "msigdbr", "clusterProfiler", "gt", "webshot2", "tidytext", "rsvg",  
          "dplyr", "VennDiagram", "stringr", "tibble", "MASS", "ggforce", "grid",
          'DiagrammeRsvg')
invisible(lapply(packs, library, character.only = TRUE))

known1 <- read.csv("step11e_PCA_main_rsIDs_known.txt", header = T, sep = "\t")
known2 <- read.csv("step11e_logistic_PCA_main_rsIDs_known.txt", header = T, sep = "\t")
known3 <- read.csv("step11e_NN_main_rsIDs_known.txt", header = T, sep = "\t")
known = union(union(known1, known2), known3)

novel1 <- read.csv("step11e_PCA_main_rsIDs_novel.txt", header = T, sep = "\t")
novel2 <- read.csv("step11e_logistic_PCA_main_rsIDs_novel.txt", header = T, sep = "\t")
novel3 <- read.csv("step11e_NN_main_rsIDs_novel.txt", header = T, sep = "\t")
novel = union(union(novel1, novel2), novel3)

#-------------------------------------------------------------------------------
# 
# NOVEL SNPS
#
#-------------------------------------------------------------------------------

my_SNPs <- read.csv("rsIDs_main_effects.txt", header = T, sep = "\t") %>%
merge(novel, by = c("rsID"))
my_SNP_gene_pairs <- read.csv("annov_main.txt", header = T, sep = "\t") %>% 
tidyr::separate("uniqID", c("chr", "pos", NA, NA)) %>% 
merge(my_SNPs, by = c("chr", "pos"))

my_gene_IDs <- my_SNP_gene_pairs[my_SNP_gene_pairs["dist"] < 300000, ]$symbol

GO_genes <- msigdbr(species = "human", subcategory = "BP")
GO_t2g <- GO_genes %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
out_GO <- enricher(gene = unique(my_gene_IDs), TERM2GENE = GO_t2g)
df_GO <- out_GO@result[out_GO@result[6] < 0.05, c(1, 6, 9, 8)]
df_GO$ID <- tolower(gsub("_", " ", df_GO$ID))

REG_genes <- msigdbr(species = "human", category = "C3")
REG_t2g <- REG_genes %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
out_REG <- enricher(gene = unique(my_gene_IDs), TERM2GENE = REG_t2g)
df_REG <- out_REG@result[out_REG@result[6] < 0.05, c(1, 6, 9, 8)]
df_REG$ID <- tolower(gsub("_", "-", df_REG$ID))

novel_table <- rbind(df_REG[c(2, 4), ], df_GO[c(1,2), ])
novel_table <- novel_table[order(novel_table$p.adjust), ]
novel_table$geneID <- gsub("/", " ", novel_table$geneID)
novel_table <- gt(novel_table, 
                  rowname_col = "ID") %>% 
  cols_label(
    p.adjust = md("Adjusted p−value"),
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
  tab_stubhead(label = "gene set ID") %>% 
  cols_width(
    ID ~ px(190),
    p.adjust ~ px(120),
    Count ~ px(100),
    geneID ~ px(170)
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
novel_table
gtsave(novel_table, "table1c.png")


#-------------------------------------------------------------------------------
# 
# KNOWN SNPS
#
#-------------------------------------------------------------------------------

my_SNPs <- read.csv("rsIDs_main_effects.txt", header = T, sep = "\t") %>%
  merge(known, by = c("rsID"))
my_SNP_gene_pairs <- read.csv("annov_main.txt", header = T, sep = "\t") %>% 
  tidyr::separate("uniqID", c("chr", "pos", NA, NA)) %>% 
  merge(my_SNPs, by = c("chr", "pos"))

my_gene_IDs <- my_SNP_gene_pairs[my_SNP_gene_pairs["dist"] < 300000, ]$symbol

GO_genes <- msigdbr(species = "human", subcategory = "BP")
GO_t2g <- GO_genes %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
out_GO <- enricher(gene = unique(my_gene_IDs), TERM2GENE = GO_t2g)
df_GO <- out_GO@result[out_GO@result[6] < 0.05, c(1, 6, 9, 8)]
df_GO$ID <- tolower(gsub("_", " ", df_GO$ID))

REG_genes <- msigdbr(species = "human", category = "C3")
REG_t2g <- REG_genes %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
out_REG <- enricher(gene = unique(my_gene_IDs), TERM2GENE = REG_t2g)
df_REG <- out_REG@result[out_REG@result[6] < 0.05, c(1, 6, 9, 8)]
df_REG$ID <- tolower(gsub("_", "-", df_REG$ID))

known_table <- df_GO[c(1,2,3,4,5,6,7,10,12,14), ]
known_table <- known_table[order(known_table$p.adjust), ]
known_table$geneID <- gsub("/", " ", known_table$geneID)
known_table <- gt(known_table, 
                  rowname_col = "ID") %>% 
  cols_label(
    p.adjust = md("Adjusted p−value"),
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
  tab_stubhead(label = "gene set ID") %>% 
  cols_width(
    ID ~ px(180),
    p.adjust ~ px(120),
    Count ~ px(95),
    geneID ~ px(275)
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
known_table
gtsave(known_table, "table1b.png")

table1a <- read.csv("table1a.txt", header = T, sep = "\t")
table1a <- gt(table1a, 
              rowname_col = "latent.phenotype.model") %>% 
  cols_label(
    Main_Novel_SNP = md("novel main effects"),
    Main_Known_SNP = md("known main effects")
  ) %>%
  cols_align(
    columns = c("Main_Novel_SNP", "Main_Known_SNP"),
    align = "right"
  ) %>% 
  tab_stubhead(label = "latent phenotype model") %>% 
  cols_width(
    latent.phenotype.model ~ px(116),
    Main_Novel_SNP ~ px(90),
    Main_Known_SNP ~ px(90)
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
table1a 
gtsave(table1a, "table1a.png")

