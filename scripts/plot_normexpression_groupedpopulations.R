### Load data

source("./scripts/load_limma.R")

# Load packages

p_load(viridis,
       palettetown,
       seriation)

list_Bioc_Pkg <- c("ComplexHeatmap",
                   "InteractiveComplexHeatmap",
                   "RColorBrewer",
                   "circlize")

BiocManager::install(pkgs = list_Bioc_Pkg,
                     update = FALSE)

## Load Bioconductor packages with pacman

invisible(lapply(list_Bioc_Pkg, function(x)
  library(x, character.only = TRUE)))

# Extract metadata and expression
## Summarise to condition level (passage, day, treatment)
## Summarise to biological replicate level, then to condition level
## Sort by:
## Passage -> Day -> cell pop -> Treatment -> biological_replicate

df_small <- df_big %>%
  group_by(
    Plate,
    cell_line,
    passage,
    day,
    treatment,
    treatment_group,
    biological_replicate,
    Detector
  ) %>%
  summarise(dCt_nexp_log2 = mean(dCt_nexp_log2, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(passage,
          day,
          cell_line,
          treatment,
          treatment_group,
          biological_replicate)

df_small <- df_small %>%
  group_by(
    passage,
    day,
    treatment_group,
    Detector
  ) %>%
  summarise(dCt_nexp_log2 = mean(dCt_nexp_log2, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(passage,
          day,
          treatment_group
          )

df_matrix <- df_small %>%
  mutate(
    col_name = str_c(
      passage,
      day,
      treatment_group,
      sep = "."
    )
  ) %>%
  dplyr::select(Detector, dCt_nexp_log2, col_name) %>%
  pivot_wider(names_from = col_name,
              values_from = dCt_nexp_log2)

matrix_exp <- as.matrix(df_matrix[-1])

row.names(matrix_exp) <- as_vector(df_matrix[1])

matrix_exp <- matrix_exp[sort(rownames(matrix_exp)), ]

df_targets <- as_tibble(colnames(matrix_exp)) %>%
  mutate(sample = value) %>%
  separate(
    value,
    into = c(
      "passage",
      "day",
      "treatment_group"
    ),
    sep = "\\."
  ) %>%
  mutate(
    sample_biol = str_c(
      passage,
      day,
      treatment_group,
      sep = "_"
    )
  ) %>%
  mutate(treatment_pair = str_c(passage, day, treatment_group, sep = "_")) %>%
  mutate(treatment_pair = factor(treatment_pair,
                                 levels = unique(treatment_pair))) %>%
  mutate(passage = factor(passage,
                          levels = c("5", "7", "13"))) %>%
  mutate(treatment_group = factor(treatment_group,
                                  levels = c("Ctrl", "50mM"))) %>%
  mutate(day = factor(day,
                      levels = c("3", "5"))) %>%
  relocate(sample)

## Create EList object

dge2 <- new("EList")
dge2$E <- matrix_exp
dge2$targets <- df_targets

genes <- tibble(name = rownames(dge2$E))
genes$group <-
  map(genes$name, \(y) names(which(sapply(gene_groups, function(x)
    y %in% x)))) %>%
  factor(levels = names(gene_groups))

genes <- genes %>%
  arrange(group, name)

dge2$genes <- genes
dge2$E     <-
  dge2$E[order(match(rownames(dge2$E), dge2$genes$name)),]

E_scaled <- t(scale(t(dge2$E)))

# Set color scheme and breaks

## For gene expression (z-scores)

col <- inferno(n = 100)
breaks <- seq(-2, 2, length.out = 100)

## For annotation data

col_treatment_group <- palettetown::pokepal(283)[c(11, 5)]

col_passage <- palettetown::pokepal(191)[c(8, 3, 9)]

col_day <- palettetown::pokepal(283)[c(3, 1)]

# Build annotation; inlclude only necessary metadata

## Dataframe of annotation data

anno <- tibble(
  `Treatment group` = dge2$targets$treatment_group,
  `Passage` = dge2$targets$passage,
  `Day` = dge2$targets$day
)

## Colour mapping

anno_cols <- list(
  `Treatment group` = c('Ctrl' = col_treatment_group[1], '50mM' = col_treatment_group[2]),
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
    `Treatment group` = list(
      nrow = length(anno_cols$`Treatment group`),
      title = 'Treatment group',
      title_position = 'topleft',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 12, fontface = 'bold'),
      labels_gp = gpar(fontsize = 12, fontface = 'bold')
    ),
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


# Create heatmap object

heatmap <- Heatmap(
  E_scaled,
  name = 'Gene\nZ-\nscore',
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
  
  # row (gene) parameters=
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

column_split  = dge2$targets$treatment_pair,
column_gap = unit(
  rep(c(0, 1),
                0.5 * length(
                   unique(dge2$targets$treatment_pair)
                   )
  ), "mm"),
  cluster_columns = F,
  show_column_dend = T,
  column_title = NULL,
  column_title_side = 'bottom',
  column_title_gp = gpar(fontsize = 12, fontface = 'bold'),
  column_title_rot = 0,
  show_column_names = FALSE,
  column_names_gp = gpar(fontsize = 10, fontface = 'bold'),
  column_names_max_height = unit(10, 'cm'),
  column_dend_height = unit(25, 'mm'),
  
  # specify top and bottom annotations
  
  top_annotation = anno_object
)

# Export heatmap of raw data

png(
  file = "./plots/heatmaps/normalised_expression_groupedpopulations.png",
  width = 12,
  height = 8,
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


# Create heatmap object - cluster everything

heatmap <- Heatmap(
  E_scaled,
  name = 'Gene\nZ-\nscore',
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
  
  # row (gene) parameters=
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
  
  # column_split  = dge2$targets$treatment_pair,
  # column_gap = unit(
  #   rep(c(0, 1),
  #       0.5 * length(
  #         unique(dge2$targets$treatment_pair)
  #       )
  #   ), "mm"),
  cluster_columns = T,
  show_column_dend = T,
  column_title = NULL,
  column_title_side = 'bottom',
  column_title_gp = gpar(fontsize = 12, fontface = 'bold'),
  column_title_rot = 0,
  show_column_names = FALSE,
  column_names_gp = gpar(fontsize = 10, fontface = 'bold'),
  column_names_max_height = unit(10, 'cm'),
  column_dend_height = unit(25, 'mm'),
  
  # specify top and bottom annotations
  
  top_annotation = anno_object
)

# Export heatmap of raw data

png(
  file = "./plots/heatmaps/normalised_expression_clustered_groupedpopulations.png",
  width = 12,
  height = 8,
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

# Create heatmap object - cluster genes only

heatmap <- Heatmap(
  E_scaled,
  name = 'Gene\nZ-\nscore',
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
  
  # row (gene) parameters=
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
  
  column_split  = dge2$targets$treatment_pair,
  column_gap = unit(
    rep(c(0, 1),
        0.5 * length(
          unique(dge2$targets$treatment_pair)
        )
    ), "mm"),
  cluster_columns = F,
  show_column_dend = T,
  column_title = NULL,
  column_title_side = 'bottom',
  column_title_gp = gpar(fontsize = 12, fontface = 'bold'),
  column_title_rot = 0,
  show_column_names = FALSE,
  column_names_gp = gpar(fontsize = 10, fontface = 'bold'),
  column_names_max_height = unit(10, 'cm'),
  column_dend_height = unit(25, 'mm'),
  
  # specify top and bottom annotations
  
  top_annotation = anno_object
)


# Export heatmap of raw data

png(
  file = "./plots/heatmaps/normalised_expression_genes_groupedpopulations.png",
  width = 12,
  height = 8,
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

