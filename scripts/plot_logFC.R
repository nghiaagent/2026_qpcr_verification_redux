### Load data

source("./scripts/analyse_limma.R")

# Load packages
# Extract metadata and logFCs

## Get gene list

fit_small <- fit_contrasts[fit_contrasts$genes$GENENAME %in% genes_GOI,]
dge_small <- dge[dge$genes$GENENAME %in% genes_GOI,]

genes <- tibble(name = rownames(dge_small$E))
genes$group <-
  map(genes$name, \(y) names(which(sapply(gene_groups, function(x)
    y %in% x)))) %>%
  factor(levels = names(gene_groups))

genes <- genes %>%
  arrange(group, name)

## Get list of logFCs

fit_small <- fit_contrasts[fit_contrasts$genes$GENENAME %in% genes_GOI,]

logFC_P5_D3 <- topTable(fit_small, coef = 1,
                        number = 100) %>%
  mutate(Detector = rownames(.)) %>%
  arrange(Detector) %>%
  mutate(passage = rep(c("5"), nrow(.))) %>%
  mutate(day     = rep(c("3"), nrow(.))) %>%
  mutate(condition = str_c(passage, day, sep = "_")) %>%
  select(condition, passage, day, Detector, logFC, adj.P.Val)

logFC_P5_D5 <- topTable(fit_small, coef = 4,
                        number = 100) %>%
  mutate(Detector = rownames(.)) %>%
  arrange(Detector) %>%
  mutate(passage = rep(c("5"), nrow(.))) %>%
  mutate(day     = rep(c("5"), nrow(.))) %>%
  mutate(condition = str_c(passage, day, sep = "_")) %>%
  select(condition, passage, day, Detector, logFC, adj.P.Val)

logFC_P7_D3 <- topTable(fit_small, coef = 2,
                        number = 100) %>%
  mutate(Detector = rownames(.)) %>%
  arrange(Detector) %>%
  mutate(passage = rep(c("7"), nrow(.))) %>%
  mutate(day     = rep(c("3"), nrow(.))) %>%
  mutate(condition = str_c(passage, day, sep = "_")) %>%
  select(condition, passage, day, Detector, logFC, adj.P.Val)

logFC_P7_D5 <- topTable(fit_small, coef = 5,
                        number = 100) %>%
  mutate(Detector = rownames(.)) %>%
  arrange(Detector) %>%
  mutate(passage = rep(c("7"), nrow(.))) %>%
  mutate(day     = rep(c("5"), nrow(.))) %>%
  mutate(condition = str_c(passage, day, sep = "_")) %>%
  select(condition, passage, day, Detector, logFC, adj.P.Val)

logFC_P13_D3 <- topTable(fit_small, coef = 3,
                         number = 100) %>%
  mutate(Detector = rownames(.)) %>%
  arrange(Detector) %>%
  mutate(passage = rep(c("13"), nrow(.))) %>%
  mutate(day     = rep(c("3"), nrow(.))) %>%
  mutate(condition = str_c(passage, day, sep = "_")) %>%
  select(condition, passage, day, Detector, logFC, adj.P.Val)

logFC_P13_D5 <- topTable(fit_small, coef = 6,
                         number = 100) %>%
  mutate(Detector = rownames(.)) %>%
  arrange(Detector) %>%
  mutate(passage = rep(c("13"), nrow(.))) %>%
  mutate(day     = rep(c("5"), nrow(.))) %>%
  mutate(condition = str_c(passage, day, sep = "_")) %>%
  select(condition, passage, day, Detector, logFC, adj.P.Val)

df_fc <- rbind(logFC_P5_D3,
               logFC_P5_D5,
               logFC_P7_D3,
               logFC_P7_D5,
               logFC_P13_D3,
               logFC_P13_D5) %>%
  mutate(passage   = factor(passage,
                            levels = c("5", "7", "13"))) %>%
  mutate(day       = factor(day,
                            levels = c("3", "5"))) %>%
  mutate(condition = factor(
    condition,
    levels = c("5_3", "5_5",
               "7_3", "7_5",
               "13_3", "13_5")
  )) %>%
  arrange(passage, day)

df_matrix_E <- df_fc %>%
  dplyr::select(Detector, logFC, condition) %>%
  pivot_wider(names_from = condition,
              values_from = logFC)

df_matrix_P <- df_fc %>%
  dplyr::select(Detector, adj.P.Val, condition) %>%
  pivot_wider(names_from = condition,
              values_from = adj.P.Val)

matrix_E <- as.matrix(df_matrix_E[-1])

matrix_P <- as.matrix(df_matrix_P[-1])

row.names(matrix_E) <- as_vector(df_matrix_E[1])

row.names(matrix_P) <- as_vector(df_matrix_P[1])

df_targets <- as_tibble(colnames(matrix_E)) %>%
  mutate(condition = value) %>%
  separate(value,
           into = c("passage", "day"),
           sep = "_") %>%
  mutate(passage   = factor(passage,
                            levels = c("5", "7", "13"))) %>%
  mutate(day       = factor(day,
                            levels = c("3", "5"))) %>%
  mutate(condition = factor(
    condition,
    levels = c("5_D3", "5_D5",
               "7_D3", "7_D5",
               "13_D3", "13_D5")
  )) %>%
  relocate(condition)

## Create EList object

dge2 <- new("EList")
dge2$E <- matrix_E
dge2$P <- matrix_P
dge2$targets <- df_targets


dge2$genes <- genes
dge2$E     <-
  dge2$E[order(match(rownames(dge2$E), dge2$genes$name)),]

dge2$P     <-
  dge2$P[order(match(rownames(dge2$P), dge2$genes$name)),]

E_scaled <- dge2$E

# Set color scheme and breaks

## For gene expression (z-scores)

#col <- turbo(n = 100)

col <- RColorBrewer::brewer.pal(10, "RdBu")

breaks <- seq(4, -4, length.out = 10)

## For annotation data

col_passage <- palettetown::pokepal(191)[c(8, 3, 9)]

col_day <- palettetown::pokepal(283)[c(3, 1)]

# Build annotation; inlclude only necessary metadata

## Dataframe of annotation data

anno <- tibble(
  `Passage` = dge2$targets$passage,
  `Day` = dge2$targets$day
)

## Colour mapping

anno_cols <- list(
  `Passage` = c('5' = col_passage[1], '7' = col_passage[2], '13' = col_passage[3]),
  `Day` = c('3' = col_day[1], '5' = col_day[2])
)

## ComplexHeatmap metadata annotation object

anno_object <- HeatmapAnnotation(
  df = anno,
  which = 'col',
  na_col = 'black',
  col = anno_cols,
  annotation_height = 0.6,
  annotation_width = unit(1, 'cm'),
  gap = unit(1, 'mm'),
  annotation_legend_param = list(
    `Passage` = list(
      nrow = length(anno_cols$`Passage`),
      title = 'Passage',
      title_position = 'topleft',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 12, fontface = 'bold'),
      labels_gp = gpar(fontsize = 12, fontface = 'bold')
    ),
    `Day` = list(
      nrow = length(anno_cols$`Day`),
      title = 'Day',
      title_position = 'topleft',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 12, fontface = 'bold'),
      labels_gp = gpar(fontsize = 12, fontface = 'bold')
    )
  )
)

# Create heatmap object - no clustering

heatmap <- Heatmap(
  E_scaled,
  name = 'logFC',
  col = colorRamp2(breaks, col),
  na_col = "grey50",
  border = T,
  
  # parameters for the colour-bar that represents gradient of expression
  
  heatmap_legend_param = list(
    color_bar = 'continuous',
    legend_direction = 'vertical',
    legend_width = unit(8, 'cm'),
    legend_height = unit(5.0, 'cm'),
    title_position = 'topcenter',
    title_gp = gpar(fontsize = 12, fontface = 'bold'),
    labels_gp = gpar(fontsize = 12, fontface = 'bold')
  ),
  
  # row (gene) parameters
  
  # row_order = dge2$genes$name,
  row_split = dge2$genes$group,
  row_gap = unit(1, "mm"),
  cluster_rows = F,
  show_row_dend = TRUE,
  cluster_row_slices = FALSE,
  #row_title = 'Statistically significant genes',
  row_title_side = 'right',
  row_title_gp = gpar(fontsize = 10,  fontface = 'bold'),
  row_title_rot = 0,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 10, fontface = 'bold'),
  row_names_side = 'right',
  row_dend_width = unit(25, 'mm'),
  
  # column (sample) parameters
  
  cluster_columns = FALSE,
  show_column_dend = FALSE,
  column_title = NULL,
  column_title_side = 'bottom',
  column_title_gp = gpar(fontsize = 12, fontface = 'bold'),
  column_title_rot = 0,
  show_column_names = FALSE,
  column_names_gp = gpar(fontsize = 10, fontface = 'bold'),
  column_names_max_height = unit(10, 'cm'),
  column_dend_height = unit(25, 'mm'),
  
  # specify top and bottom annotations
  
  top_annotation = anno_object,
  
  # Add p-value to cells
  cell_fun = function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%.3f", dge2$P[i, j]), x, y, gp = gpar(fontsize = 10))
  }
)


# Export heatmap of raw data

png(
  file = "./plots/heatmaps/logFC_Pval.png",
  width = 12,
  height = 6,
  units = "in",
  res = 1200
)

heatmap <- draw(
  heatmap,
  heatmap_legend_side = 'left',
  annotation_legend_side = 'right',
  row_sub_title_side = 'left'
)

dev.off()


# Create heatmap object - cluster genes, group kept

heatmap <- Heatmap(
  E_scaled,
  name = 'logFC',
  col = colorRamp2(breaks, col),
  na_col = "grey50",
  border = T,
  
  # parameters for the colour-bar that represents gradient of expression
  
  heatmap_legend_param = list(
    color_bar = 'continuous',
    legend_direction = 'vertical',
    legend_width = unit(8, 'cm'),
    legend_height = unit(5.0, 'cm'),
    title_position = 'topcenter',
    title_gp = gpar(fontsize = 12, fontface = 'bold'),
    labels_gp = gpar(fontsize = 12, fontface = 'bold')
  ),
  
  # row (gene) parameters
  
  # row_order = dge2$genes$name,
  row_split = dge2$genes$group,
  row_gap = unit(1, "mm"),
  cluster_rows = T,
  show_row_dend = TRUE,
  cluster_row_slices = FALSE,
  #row_title = 'Statistically significant genes',
  row_title_side = 'right',
  row_title_gp = gpar(fontsize = 10,  fontface = 'bold'),
  row_title_rot = 0,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 10, fontface = 'bold'),
  row_names_side = 'right',
  row_dend_width = unit(25, 'mm'),
  
  # column (sample) parameters
  
  cluster_columns = FALSE,
  show_column_dend = FALSE,
  column_title = NULL,
  column_title_side = 'bottom',
  column_title_gp = gpar(fontsize = 12, fontface = 'bold'),
  column_title_rot = 0,
  show_column_names = FALSE,
  column_names_gp = gpar(fontsize = 10, fontface = 'bold'),
  column_names_max_height = unit(10, 'cm'),
  column_dend_height = unit(25, 'mm'),
  
  # specify top and bottom annotations
  
  top_annotation = anno_object,
  
  # Add p-value to cells
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.3f", dge2$P[i, j]), x, y, gp = gpar(fontsize = 10))
  }
)


# Export heatmap of raw data

png(
  file = "./plots/heatmaps/logFC_clustered_Pval.png",
  width = 12,
  height = 6,
  units = "in",
  res = 1200
)

heatmap <- draw(
  heatmap,
  heatmap_legend_side = 'left',
  annotation_legend_side = 'right',
  row_sub_title_side = 'left'
)

dev.off()

# Create heatmap object - no clustering

heatmap <- Heatmap(
  E_scaled,
  name = 'logFC',
  col = colorRamp2(breaks, col),
  na_col = "grey50",
  border = T,
  
  # parameters for the colour-bar that represents gradient of expression
  
  heatmap_legend_param = list(
    color_bar = 'continuous',
    legend_direction = 'vertical',
    legend_width = unit(8, 'cm'),
    legend_height = unit(5.0, 'cm'),
    title_position = 'topcenter',
    title_gp = gpar(fontsize = 12, fontface = 'bold'),
    labels_gp = gpar(fontsize = 12, fontface = 'bold')
  ),
  
  # row (gene) parameters
  
  # row_order = dge2$genes$name,
  # row_split = dge2$genes$group,
  row_gap = unit(1, "mm"),
  cluster_rows = T,
  show_row_dend = TRUE,
  cluster_row_slices = FALSE,
  #row_title = 'Statistically significant genes',
  row_title_side = 'right',
  row_title_gp = gpar(fontsize = 10,  fontface = 'bold'),
  row_title_rot = 0,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 10, fontface = 'bold'),
  row_names_side = 'right',
  row_dend_width = unit(25, 'mm'),
  
  # column (sample) parameters
  
  cluster_columns = FALSE,
  show_column_dend = FALSE,
  column_title = NULL,
  column_title_side = 'bottom',
  column_title_gp = gpar(fontsize = 12, fontface = 'bold'),
  column_title_rot = 0,
  show_column_names = FALSE,
  column_names_gp = gpar(fontsize = 10, fontface = 'bold'),
  column_names_max_height = unit(10, 'cm'),
  column_dend_height = unit(25, 'mm'),
  
  # specify top and bottom annotations
  
  top_annotation = anno_object,
  
  # Add p-value to cells
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.3f", dge2$P[i, j]), x, y, gp = gpar(fontsize = 10))
  }
)


# Export heatmap of raw data

png(
  file = "./plots/heatmaps/logFC_clustered_nogroup_Pval.png",
  width = 12,
  height = 6,
  units = "in",
  res = 1200
)

heatmap <- draw(
  heatmap,
  heatmap_legend_side = 'left',
  annotation_legend_side = 'right',
  row_sub_title_side = 'left'
)

dev.off()

