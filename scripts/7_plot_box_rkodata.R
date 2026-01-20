# Import packages
library("DESeq2")
library("tidyverse")

# Set up axes label mapping
labels_map_rkodata <- c(
  "5 3 Ctrl" = "P5 D3 Ctrl",
  "5 3 10ug" = "P5 D3 Hep",
  "5 5 Ctrl" = "P5 D5 Ctrl",
  "5 5 10ug" = "P5 D5 Hep",
  "7 3 Ctrl" = "P7 D3 Ctrl",
  "7 3 10ug" = "P7 D3 Hep",
  "7 5 Ctrl" = "P7 D5 Ctrl",
  "7 5 10ug" = "P7 D5 Hep",
  "13 3 Ctrl" = "P13 D3 Ctrl",
  "13 3 10ug" = "P13 D3 Hep",
  "13 5 Ctrl" = "P13 D5 Ctrl",
  "13 5 10ug" = "P13 D5 Hep"
)
labels_map_txome <- c(
  "P5D3Untreated" = "P5 D3 Ctrl",
  "P5D3Treated" = "P5 D3 Hep",
  "P5D5Untreated" = "P5 D5 Ctrl",
  "P5D5Treated" = "P5 D5 Hep",
  "P7D3Untreated" = "P7 D3 Ctrl",
  "P7D3Treated" = "P7 D3 Hep",
  "P7D5Untreated" = "P7 D5 Ctrl",
  "P7D5Treated" = "P7 D5 Hep",
  "P13D3Untreated" = "P13 D3 Ctrl",
  "P13D3Treated" = "P13 D3 Hep",
  "P13D5Untreated" = "P13 D5 Ctrl",
  "P13D5Treated" = "P13 D5 Hep"
)

# Set genes for exclusion and GOIs
## RKO data
rkodata_genes_interest <- c(
  "GPC1",
  "GPC4",
  "SDC1",
  "SDC2",
  "SDC3",
  "SDC4",
  "EXT1",
  "EXT2",
  "NDST1",
  "NDST2",
  "HS2ST1",
  "ACTA2",
  "CD44",
  "NES",
  "OCT3/4",
  "TUBB3",
  "ENO2"
)

rkodata_genes_exclude <- c(
  "ALB",
  "COL2A1",
  "FGFR2",
  "FOXA2",
  "GFAP",
  "GPC2",
  "GPC3",
  "GPC5",
  "GPC6",
  "HPSE",
  "HS6ST2",
  "HS6ST3",
  "MAP2",
  "MSI1",
  "NEFM",
  "NDST3",
  "NDST4",
  "NEUROG2",
  "SOX1",
  "RUNX2",
  "VTN"
)

## Txome data
txome_genes_interest <- c(
  "GPC1",
  "GPC4",
  "SDC1",
  "SDC2",
  "SDC3",
  "SDC4",
  "EXT1",
  "EXT2",
  "NDST1",
  "NDST2",
  "HS2ST1",
  "ACTA2",
  "CD44",
  "NES",
  "POU5F1",
  "TUBB3",
  "ENO2"
)

# Load data
## RKO data
### df
## Filter for only experiment of interest
## Select only necessary columns
## Add columns for condition
## Rename to match current data format
## Recode to factor where appropriate
rkodata_df_big <- readRDS(
  here::here(
    "input",
    "rkodata",
    "df_big.rds"
  )
)

rkodata_df_big <- rkodata_df_big %>%
  dplyr::filter(
    !Detector %in% rkodata_genes_exclude,
    cell_line %in% c("hMSC-20176", "hMSC-21558"),
    day %in% c("3"),
    treatment == "Hep"
  ) %>%
  dplyr::mutate(
    condition_ID = str_c(
      passage,
      day,
      treatment_group,
      sep = " "
    )
  ) %>%
  dplyr::select(
    sample,
    Plate,
    condition_ID,
    biological_replicate,
    technical_replicate,
    cell_line,
    passage,
    day,
    treatment_group,
    Detector,
    Ct,
    Ct_18S,
    dCt,
    dCt_neg,
    dCt_nexp_log2,
    dCt_nexp_log2_offset
  ) %>%
  dplyr::rename(
    "name" = "sample",
    "filename" = "Plate",
    "condition_id" = "condition_ID",
    "treatment" = "treatment_group",
    "gene" = "Detector",
    "ct" = "Ct",
    "ct_18S" = "Ct_18S",
    "dct" = "dCt",
    "dct_neg" = "dCt_neg",
    "dct_nexp_log2" = "dCt_nexp_log2",
    "dct_nexp_log2_offset" = "dCt_nexp_log2_offset"
  ) %>%
  dplyr::arrange(
    cell_line,
    passage,
    day,
    treatment,
    gene,
    biological_replicate,
    technical_replicate
  ) %>%
  dplyr::mutate(
    condition_id = condition_id %>%
      factor() %>%
      fct_inorder()
  ) %>%
  droplevels()

### Limma fit of heparin data
rkodata_fit_contrasts <- readRDS(
  here::here(
    "input",
    "rkodata",
    "4_fit_contrasts.rds"
  )
)

### Filter for experiment of interest
rkodata_fit_contrasts <- rkodata_fit_contrasts[["heparin"]]

## txome data
txome_quant_deseq2 <- readRDS(
  file = here::here(
    "input",
    "txome",
    "quant_deseq2_batchcor_nofilter.RDS"
  )
)

txome_results <- readRDS(
  file = here::here(
    "input",
    "txome",
    "results_deseq2_nofilter.RDS"
  )
)
### Subset data based on conditions and genes
txome_quant_small <- txome_quant_deseq2 %$%
  .[,
    colData(.)$condition_ID %in%
      c(
        "P5D3Untreated",
        "P5D3Treated",
        "P7D3Untreated",
        "P7D3Treated",
        "P13D3Untreated",
        "P13D3Treated"
      )
  ]

txome_genes_interest_id <- txome_quant_small %>%
  rowRanges() %>%
  as.data.frame() %>%
  filter(symbol %in% txome_genes_interest) %>%
  arrange(match(gene_name, txome_genes_interest))

### Readjust DESeq results to include only GOIs
results_subset <- txome_results %>%
  map(\(x) {
    x <- x[rownames(x) %in% txome_genes_interest_id$gene_id, ]
    x$padj <- p.adjust(x$pvalue, method = "BH")
    return(x)
  })

# Define functions for drawing boxplots
## RKO data
rkodata_plot_goi <- function(gene_sel) {
  # Get table containing dCt for hMSC
  quant_small <- rkodata_df_big %>%
    dplyr::filter(gene %in% gene_sel) %>%
    droplevels()

  # Get position to start drawing signif bars
  y_position <- max(quant_small$dct_nexp_log2) * 0.95

  # Plot data
  plot <- quant_small %>%
    # Define plot x y axes
    ggplot(aes(x = condition_id, y = dct_nexp_log2)) +
    geom_boxplot(aes(color = condition_id)) +
    geom_jitter(aes(color = condition_id)) +
    scale_x_discrete(labels = labels_map_rkodata) +
    scale_y_continuous(expand = expansion(0, 0.3)) +
    scale_colour_manual(
      values = palette_merge[c(1:6)],
      labels = c(
        "5 3 Ctrl" = "P5 D3 Ctrl",
        "5 3 10ug" = "P5 D3 Hep",
        "7 3 Ctrl" = "P7 D3 Ctrl",
        "7 3 10ug" = "P7 D3 Hep",
        "13 3 Ctrl" = "P13 D3 Ctrl",
        "13 3 10ug" = "P13 D3 Hep"
      )
    ) +
    theme_classic() +
    ggtitle(
      label = gene_sel %>%
        str_replace_all("_", " ") %>%
        str_wrap(
          width = 20,
          whitespace_only = TRUE
        )
    ) +
    ylab("Norm. log expression")

  # Return object
  return(plot)
}

## TXome data
txome_plot_goi <- function(gene_sel, gene_name) {
  # Get table containing normalised gene counts
  gene_counts <- txome_quant_small %>%
    plotCounts(
      gene = gene_sel,
      intgroup = c("condition_ID"),
      returnData = TRUE
    ) %>%
    droplevels()

  # Get position to start drawing signif bars
  y_position <- log10(max(gene_counts$count) * 1.1)

  # Plot data
  plot <- gene_counts %>%
    ggplot(aes(x = condition_ID, y = count)) +
    geom_boxplot(aes(color = condition_ID)) +
    geom_jitter(aes(color = condition_ID)) +
    scale_x_discrete(labels = labels_map_txome) +
    scale_y_log10(expand = expansion(0, 0.05)) +
    scale_colour_manual(
      values = palette_merge[c(1, 2, 3, 4, 5, 6)],
      labels = c(
        "P5D3Untreated" = "P5 D3 Ctrl",
        "P5D3Treated" = "P5 D3 Hep",
        "P7D3Untreated" = "P7 D3 Ctrl",
        "P7D3Treated" = "P7 D3 Hep",
        "P13D3Untreated" = "P13 D3 Ctrl",
        "P13D3Treated" = "P13 D3 Hep"
      )
    ) +
    theme_classic() +
    ggtitle(label = gene_name) +
    ylab("Norm. counts")

  # Return object
  return(plot)
}

# Draw list of boxplots
## RKO data
rkodata_plots_goi <- rkodata_genes_interest %>%
  map(\(gene) rkodata_plot_goi(gene)) %>%
  set_names(rkodata_genes_interest)

## Txome data
txome_plots_goi <- map2(
  txome_genes_interest_id$gene_id,
  txome_genes_interest_id$gene_name,
  \(gene_sel, gene_name) txome_plot_goi(gene_sel, gene_name)
) %>%
  set_names(txome_genes_interest_id$gene_name)

# Save data
saveRDS(
  rkodata_plots_goi,
  here::here(
    "output",
    "plots_box",
    "data",
    "rkodata_plots_goi.rds"
  )
)

saveRDS(
  txome_plots_goi,
  here::here(
    "output",
    "plots_box",
    "data",
    "txome_plots_goi.rds"
  )
)
