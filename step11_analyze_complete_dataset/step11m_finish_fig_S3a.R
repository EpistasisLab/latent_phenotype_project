figure_S3a <- read.csv("figureS3a.txt", header = T, sep = "\t")
colnames(figure_S3a) <- gsub("\\.", " ", colnames(figure_S3a))
figure_S3a <- gt(figure_S3a, 
                 rowname_col = "ICD10 code") %>% 
  fmt_integer(
    columns = c("PCA average importance rank", 
                "Autoencoder average importance rank", 
                "Logistic PCA average importance rank")
  ) %>% 
  cols_width(
    `PCA average importance rank` ~ px(120),
    `Autoencoder average importance rank` ~ px(120),
    `Logistic PCA average importance rank` ~ px(120),
  ) %>%
  tab_stubhead(label = "ICD10 code") %>% 
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
  opt_horizontal_padding(scale = 3) %>% 
  opt_vertical_padding(scale = 0)
figure_S3a 
gtsave(figure_S3a , "figure_S3a.png")