table1a_df <- read.csv("table1a.txt", header = T, sep = "\t", check.names = FALSE)
table1b_df <- read.csv("table1b.txt", header = T, sep = "\t", check.names = FALSE)
table1c_df <- read.csv("table1c.txt", header = T, sep = "\t", check.names = FALSE)
table1d_df <- read.csv("table1d.txt", header = T, sep = "\t", check.names = FALSE)

figure3b_df <- read.csv("figure3b.txt", header = F, sep = "\t", check.names = FALSE)
colnames(figure3b_df) <- c( "rsID", "binary AHF p-value", "Aragam, et al p-value*")
figure3b_df["Aragam, et al p-value*"] <- c(9.08e-10, 2.15e-7, 2.63e-7, 2.65e-7, 3.55e-7, 1e-3, NA)

figure3b_PCA <- read.csv("figure3b_TRACE_PCA.txt", header = F, sep = "\t", check.names = FALSE)
figure3b_PCA["V9"] <- "PCA"
figure3b_NN <- read.csv("figure3b_TRACE_NN.txt", header = F, sep = "\t", check.names = FALSE)
figure3b_NN["V9"] <- "autoencoder"
figure3b_logistic_PCA <- read.csv("figure3b_TRACE_logistic_PCA.txt", header = F, sep = "\t", check.names = FALSE)
figure3b_logistic_PCA["V9"] <- "logistic PCA"

figure3b_other <- rbind(figure3b_PCA, figure3b_NN, figure3b_logistic_PCA)
figure3b_other <- figure3b_other %>%
  group_by(V1) %>%
  summarise(`best latent phenotype p-value` = min(V4, na.rm = TRUE), 
            `latent phenotype model` = V9[which.min(V4)]) %>%
  rename(rsID = V1)

figure3b_df <- merge(figure3b_other, figure3b_df, by = "rsID")
figure3b_df <- figure3b_df[order(figure3b_df[,4]), ]

table1a <- gt(table1a_df ,
              rowname_col = "latent phenotype model") %>%
  tab_options(
    column_labels.background.color = "#0b1d78",
    table.background.color = "#d8f9ff"
  ) %>%
  fmt_scientific(
    columns = c("p-value", "conditional p-value"),
    decimals = 2
  ) %>%
  fmt_number(
    columns = c("TRACE vs. linear E[logp]", "conditional E[logp]"),
    decimals = 2
  ) %>%
  cols_width(
    everything() ~ px(115)
  ) %>%
  cols_width(
    `TRACE vs. linear E[logp]` ~ px(90),
    `conditional E[logp]` ~ px(90),
  ) %>%
  cols_width(
    `latent phenotype model` ~ px(120)
  ) %>%
  tab_stubhead(label = "latent phenotype model") %>%
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
  cols_label(
    `TRACE vs. linear E[logp]` = "TRACE vs linear ΔE[logp]",
    `conditional E[logp]` = "conditional ΔE[logp]"
  ) %>%
  opt_horizontal_padding(scale = 2)
table1a
gtsave(table1a, "table1a.png")


table1b <- gt(table1b_df ,
              rowname_col = "latent phenotype model") %>%
  tab_options(
    column_labels.background.color = "#0b1d78",
    table.background.color = "#d8f9ff"
  ) %>%
  cols_width(
    everything() ~ px(75)
  ) %>%
  cols_width(
    `latent phenotype model` ~ px(120)
  ) %>%
  tab_stubhead(label = "latent phenotype model") %>%
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
  cols_label(
    GxAlcohol_Novel_SNP = "novel SNPs",
    GxAlcohol_Known_SNP = "known SNPs",
    GxGender_Novel_SNP = "novel SNPs",
    GxGender_Known_SNP = "known SNPs",
    GxSmoking_Novel_SNP = "novel SNPs",
    GxSmoking_Known_SNP = "known SNPs",
    Main_Novel_SNP = "novel SNPs",
    Main_Known_SNP = "known SNPs"
  ) %>%
  tab_spanner(
    label = "main effects",
    columns = vars(Main_Novel_SNP, Main_Known_SNP)
  ) %>%
  tab_spanner(
    label = "GxSmoking effects",
    columns = vars(GxSmoking_Novel_SNP, GxSmoking_Known_SNP)
  ) %>%
  tab_spanner(
    label = "GxGender effects",
    columns = vars(GxGender_Novel_SNP, GxGender_Known_SNP)
  ) %>%
  tab_spanner(
    label = "GxAlcohol effects",
    columns = vars(GxAlcohol_Novel_SNP, GxAlcohol_Known_SNP)
  ) %>%
opt_horizontal_padding(scale = 2)
table1b
gtsave(table1b, "table1b.png")


table1c <- gt(table1c_df ,
              rowname_col = "main model") %>%
  tab_options(
    column_labels.background.color = "#0b1d78",
    table.background.color = "#d8f9ff"
  ) %>%
  fmt_scientific(
    columns = "p_diff alt minus null",
    decimals = 2
  ) %>%
  fmt_number(
    columns = c("main E[-log10(p)]", "alt E[-log10(p)]", "null E[-log10(p)]"),
    decimals = 2
  ) %>%
  cols_width(
    everything() ~ px(120)
  ) %>%
  cols_width(
    `main E[-log10(p)]` ~ px(87),
    `alt E[-log10(p)]` ~ px(87),
    `null E[-log10(p)]` ~ px(87)
  ) %>%
  tab_stubhead(label = "main") %>%
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
  cols_label(
    `main E[-log10(p)]` = "main E[−logp]",
    `alt E[-log10(p)]` = "alt E[−logp]",
    `null E[-log10(p)]` = "null E[−logp]",
  ) %>%
  opt_horizontal_padding(scale = 2)
table1c
gtsave(table1c, "table1c.png")

table1d <- gt(table1d_df ,
              rowname_col = "main model") %>%
  tab_options(
    column_labels.background.color = "#0b1d78",
    table.background.color = "#d8f9ff"
  ) %>%
  fmt_scientific(
    columns = "p_diff alt minus null",
    decimals = 2
  ) %>%
  fmt_number(
    columns = c("main E[-log10(p)]", "alt E[-log10(p)]", "null E[-log10(p)]"),
    decimals = 2
  ) %>%
  cols_width(
    everything() ~ px(120)
  ) %>%
  cols_width(
    `main E[-log10(p)]` ~ px(87),
    `alt E[-log10(p)]` ~ px(87),
    `null E[-log10(p)]` ~ px(87)
  ) %>%
  tab_stubhead(label = "main") %>%
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
  cols_label(
    `main E[-log10(p)]` = "main E[−logp]",
    `alt E[-log10(p)]` = "alt E[−logp]",
    `null E[-log10(p)]` = "null E[−logp]",
  ) %>%
  opt_horizontal_padding(scale = 2)
table1d
gtsave(table1d, "table1d.png")

figure3b <- gt(figure3b_df ,
              rowname_col = "rsID") %>%
  tab_options(
    column_labels.background.color = "#0b1d78",
    table.background.color = "#d8f9ff"
  ) %>%
  fmt_scientific(
    columns = c("best latent phenotype p-value", 
                "binary AHF p-value",
                "Aragam, et al p-value*"),
    decimals = 2
  ) %>%
  cols_width(
    everything() ~ px(110)
  ) %>%
  tab_stubhead(label = "rsID") %>%
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
  opt_horizontal_padding(scale = 2)
figure3b
gtsave(figure3b, "figure3b.png")

# Set parameters for multivariate data
n <- 250 # number of observations
mu <- c(0, 0) # means
Sigma <- matrix(c(1, 0.5, 0.5, 0.7), 2) # covariance matrix

# Generate multivariate data
set.seed(123) # for reproducible results
data <- MASS::mvrnorm(n, mu, Sigma)
df <- data.frame(data)

# Calculate principal components
pc <- stats::princomp(df)

# Extract the PCA rotation vectors (loadings)
rotation_vectors <- 3*(pc$sdev**2)*pc$loadings[,1:2]

# Create data frame for rotation vectors
df_pc <- data.frame(rotation_vectors)
colnames(df_pc) <- c("PC1", "PC2")

# Plot the data
ggplot(df, aes(x = X1, y = X2)) +
  geom_point(aes(color = 'Observed\nPhenotype\nVectors')) +
  scale_color_manual(name = "", values = c('Observed\nPhenotype\nVectors' = 'black')) +
  geom_segment(data = df_pc, aes(x = 0, y = 0, xend = PC1, yend = PC2, linetype = 'Latent\nPhenotype\nDimensions'), 
               arrow = arrow(length = unit(0.30, "cm")), size = 1.2, color = "red") +
  scale_linetype_manual(name = "", values = c('Latent\nPhenotype\nDimensions' = "solid")) +
  coord_fixed() +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.key = element_blank(),
    legend.position = c(.24, .70),
    legend.spacing.y = unit(-0.2, "cm"),
    legend.justification = c("right", "bottom"),
    legend.box.just = "left",
    text = element_text(size = 14)
  ) +
  labs(
    title = "Relationship Between Observed\nand Latent Phenotypes",
    x = "Observable phenotype 1",
    y = "Observable phenotype 2"
  )
png(filename="figure1b.png", width=2700, height=1800, units="px", res=300)
last_plot()
dev.off()


# create the graph
p = 1.4
a = 0.7
graph <- create_graph() %>%
  add_node(label = "E", type = "E", node_aes = node_aes(shape = "circle", fill = "blue", color = "black")) %>%
  add_node(label = "G", type = "G", node_aes = node_aes(shape = "circle", fill = "blue", color = "black")) %>%
  add_node(label = "G", type = "G", node_aes = node_aes(shape = "circle", fill = "blue", color = "black")) %>%
  add_node(label = "G", type = "G", node_aes = node_aes(shape = "circle", fill = "blue", color = "black")) %>%
  add_node(label = "T", type = "T", node_aes = node_aes(shape = "circle", fill = "red", color = "black")) %>%
  add_node(label = "T", type = "T", node_aes = node_aes(shape = "circle", fill = "red", color = "black")) %>%
  add_node(label = "P", type = "P", node_aes = node_aes(shape = "circle", fill = "red", color = "black")) %>%
  add_node(label = "P", type = "P", node_aes = node_aes(shape = "circle", fill = "red", color = "black")) %>%
  add_node(label = "Z", type = "Z", node_aes = node_aes(shape = "circle", fill = "#3cb371", color = "black")) %>%
  add_node(label = "Z", type = "Z", node_aes = node_aes(shape = "circle", fill = "#3cb371", color = "black")) %>%
  add_node(label = "Y", type = "Y", node_aes = node_aes(shape = "circle", fill = "blue", color = "black")) %>%
  add_node(label = "Y", type = "Y", node_aes = node_aes(shape = "circle", fill = "blue", color = "black")) %>%
  add_node(label = "Y", type = "Y", node_aes = node_aes(shape = "circle", fill = "blue", color = "black")) %>%
  add_edge(from = 1, to = 5, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>%
  add_edge(from = 1, to = 8, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>%
  add_edge(from = 2, to = 7, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>%
  add_edge(from = 3, to = 5, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>% 
  add_edge(from = 3, to = 6, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>% 
  add_edge(from = 4, to = 6, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>% 
  add_edge(from = 5, to = 7, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>%
  add_edge(from = 6, to = 8, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>% 
  add_edge(from = 7, to = 8, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>% 
  add_edge(from = 7, to = 9, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>% 
  add_edge(from = 8, to = 10, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>% 
  add_edge(from = 9, to = 11, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>% 
  add_edge(from = 10, to = 11, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>% 
  add_edge(from = 9, to = 12, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>% 
  add_edge(from = 10, to = 12, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>% 
  add_edge(from = 12, to = 13, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>% 
  set_node_position(node = 1, x = 2, y = 1) %>% 
  set_node_position(node = 2, x = 4, y = 2) %>% 
  set_node_position(node = 3, x = 2, y = 2) %>% 
  set_node_position(node = 4, x = 0, y = 2) %>% 
  set_node_position(node = 5, x = 3, y = 3) %>% 
  set_node_position(node = 6, x = 1, y = 3) %>% 
  set_node_position(node = 7, x = 3, y = 4) %>% 
  set_node_position(node = 8, x = 1, y = 4) %>% 
  set_node_position(node = 9, x = 3, y = 5) %>% 
  set_node_position(node = 10, x = 1, y = 5) %>% 
  set_node_position(node = 11, x = 3, y = 6) %>% 
  set_node_position(node = 12, x = 1, y = 6) %>% 
  set_node_position(node = 13, x = 0, y = 6)

render_graph(graph)

p = 1.4
a = 0.7
graph <- create_graph() %>%
  add_node(label = "G", type = "G", node_aes = node_aes(shape = "circle", fill = "blue", color = "black")) %>%
  add_node(label = "G", type = "G", node_aes = node_aes(shape = "circle", fill = "blue", color = "black")) %>%
  add_node(label = "G", type = "G", node_aes = node_aes(shape = "circle", fill = "blue", color = "black")) %>%
  add_node(label = "G", type = "G", node_aes = node_aes(shape = "circle", fill = "blue", color = "black")) %>%
  add_node(label = "G", type = "G", node_aes = node_aes(shape = "circle", fill = "blue", color = "black")) %>%
  add_node(label = "Z", type = "Z", node_aes = node_aes(shape = "circle", fill = "#3cb371", color = "black")) %>%
  add_node(label = "Y", type = "Y", node_aes = node_aes(shape = "circle", fill = "blue", color = "black")) %>%
  add_node(label = "Y", type = "Y", node_aes = node_aes(shape = "circle", fill = "blue", color = "black")) %>%
  add_node(label = "Y", type = "Y", node_aes = node_aes(shape = "circle", fill = "blue", color = "black")) %>%
  add_node(label = "Y", type = "Y", node_aes = node_aes(shape = "circle", fill = "blue", color = "black")) %>%
  add_node(label = "Y", type = "Y", node_aes = node_aes(shape = "circle", fill = "blue", color = "black")) %>%
  add_edge(from = 1, to = 6, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>%
  add_edge(from = 2, to = 6, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>%
  add_edge(from = 3, to = 6, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>%
  add_edge(from = 4, to = 6, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>% 
  add_edge(from = 5, to = 6, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>%
  add_edge(from = 6, to = 7, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>%
  add_edge(from = 6, to = 8, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>%
  add_edge(from = 6, to = 9, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>%
  add_edge(from = 6, to = 10, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>% 
  add_edge(from = 6, to = 11, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>% 
  set_node_position(node = 1, x = 0, y = 0) %>% 
  set_node_position(node = 2, x = 1, y = 0) %>% 
  set_node_position(node = 3, x = 2, y = 0) %>% 
  set_node_position(node = 4, x = 3, y = 0) %>% 
  set_node_position(node = 5, x = 4, y = 0) %>% 
  set_node_position(node = 6, x = 2, y = 1) %>% 
  set_node_position(node = 7, x = 0, y = 2) %>% 
  set_node_position(node = 8, x = 1, y = 2) %>% 
  set_node_position(node = 9, x = 2, y = 2) %>% 
  set_node_position(node = 10, x = 3, y = 2) %>% 
  set_node_position(node = 11, x = 4, y = 2)

render_graph(graph)
# must be manually saved as figure 1a1

p = 1.4
a = 0.7
graph <- create_graph() %>%
  add_node(label = "G", type = "G", node_aes = node_aes(shape = "circle", fill = "blue", color = "black")) %>%
  add_node(label = "G", type = "G", node_aes = node_aes(shape = "circle", fill = "blue", color = "black")) %>%
  add_node(label = "G", type = "G", node_aes = node_aes(shape = "circle", fill = "blue", color = "black")) %>%
  add_node(label = "G", type = "G", node_aes = node_aes(shape = "circle", fill = "blue", color = "black")) %>%
  add_node(label = "G", type = "G", node_aes = node_aes(shape = "circle", fill = "blue", color = "black")) %>%
  add_node(label = "Y", type = "Y", node_aes = node_aes(shape = "circle", fill = "blue", color = "black")) %>%
  add_node(label = "Y", type = "Y", node_aes = node_aes(shape = "circle", fill = "blue", color = "black")) %>%
  add_node(label = "Y", type = "Y", node_aes = node_aes(shape = "circle", fill = "blue", color = "black")) %>%
  add_node(label = "Y", type = "Y", node_aes = node_aes(shape = "circle", fill = "blue", color = "black")) %>%
  add_node(label = "Y", type = "Y", node_aes = node_aes(shape = "circle", fill = "blue", color = "black")) %>%
  add_edge(from = 1, to = 6, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>%
  add_edge(from = 1, to = 7, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>%
  add_edge(from = 1, to = 8, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>%
  add_edge(from = 1, to = 9, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>% 
  add_edge(from = 1, to = 10, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>%
  add_edge(from = 2, to = 6, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>%
  add_edge(from = 2, to = 7, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>%
  add_edge(from = 2, to = 8, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>%
  add_edge(from = 2, to = 9, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>% 
  add_edge(from = 2, to = 10, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>%
  add_edge(from = 3, to = 6, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>%
  add_edge(from = 3, to = 7, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>%
  add_edge(from = 3, to = 8, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>%
  add_edge(from = 3, to = 9, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>% 
  add_edge(from = 3, to = 10, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>%
  add_edge(from = 4, to = 6, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>%
  add_edge(from = 4, to = 7, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>%
  add_edge(from = 4, to = 8, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>%
  add_edge(from = 4, to = 9, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>% 
  add_edge(from = 4, to = 10, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>%
  add_edge(from = 5, to = 6, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>%
  add_edge(from = 5, to = 7, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>%
  add_edge(from = 5, to = 8, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>%
  add_edge(from = 5, to = 9, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>% 
  add_edge(from = 5, to = 10, edge_aes = edge_aes(color = "black",  penwidth = p, arrowsize = a)) %>%
  set_node_position(node = 1, x = 0, y = 0) %>% 
  set_node_position(node = 2, x = 1, y = 0) %>% 
  set_node_position(node = 3, x = 2, y = 0) %>% 
  set_node_position(node = 4, x = 3, y = 0) %>% 
  set_node_position(node = 5, x = 4, y = 0) %>% 
  set_node_position(node = 6, x = 0, y = 1) %>% 
  set_node_position(node = 7, x = 1, y = 1) %>% 
  set_node_position(node = 8, x = 2, y = 1) %>% 
  set_node_position(node = 9, x = 3, y = 1) %>% 
  set_node_position(node = 10, x = 4, y = 1)

render_graph(graph)
# must be manually saved as figure 1a2