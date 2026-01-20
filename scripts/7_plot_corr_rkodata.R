# Import packages
library("data.table")
library("DESeq2")
library("tidyverse")

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

txome_conditions_interest <- c(
  "P5D3Untreated",
  "P5D3Treated",
  "P7D3Untreated",
  "P7D3Treated",
  "P13D3Untreated",
  "P13D3Treated"
)

# Set mapping between RKO data and txome sample names
mapping <- data.frame(
  names_rkodata = str_c(
    # Cell pops
    rep(c("hMSC-20176", "hMSC-21558"), c(36, 36)),
    # Passages
    rep(c("5", "7", "13"), c(12, 12, 12)) %>%
      rep(2),
    # Day
    rep(c("3", "5"), c(6, 6)) %>%
      rep(6),
    # Treatment
    rep(c("Ctrl", "10ug"), c(3, 3)) %>%
      rep(12),
    # Biol rep
    rep(c("1", "2", "3"), 24),
    sep = "_"
  ),
  names_txome = str_c(
    # Cell pops
    rep(c("hMSC-20176", "hMSC-21558"), c(36, 36)),
    # Passages
    rep(c("P5", "P7", "P13"), c(12, 12, 12)) %>%
      rep(2),
    # Day
    rep(c("D3", "D5"), c(6, 6)) %>%
      rep(6),
    # Treatment
    rep(c("0", "10"), c(3, 3)) %>%
      rep(12),
    sep = "_"
  ) %>%
    str_c(
      # Biol rep
      rep(c("1", "2", "3"), 24),
      sep = "-"
    )
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

## txome data
txome_quant_deseq2 <- readRDS(
  file = here::here(
    "input",
    "txome",
    "quant_deseq2_batchcor_nofilter.RDS"
  )
)

# Clean data
# Get table of per biological replicate gene expression

## qPCR data
rkodata_quant_sel <- rkodata_df_big %>%
  # Subset data based on genes and condition
  dplyr::filter(gene %in% rkodata_genes_interest) %>%
  # Summarise to biol replicate level
  dplyr::group_by(name, gene) %>%
  dplyr::summarise(
    condition_id = condition_id %>%
      unique(),
    biological_replicate = biological_replicate %>%
      unique(),
    passage = passage %>%
      unique(),
    day = day %>%
      unique(),
    treatment = treatment %>%
      unique(),
    dct_neg_mean = dct_neg %>%
      mean(na.rm = TRUE)
  ) %>%
  dplyr::mutate(
    name = name %>%
      str_replace_all(" ", "_")
  ) %>%
  droplevels()

## txome data
### Subset data based on conditions and genes
txome_quant_small <- txome_quant_deseq2 %$%
  .[, colData(.)$condition_ID %in% txome_conditions_interest]

### Get gene ids
txome_genes_sel <- txome_quant_small %>%
  rowRanges() %>%
  as.data.frame() %>%
  filter(symbol %in% txome_genes_interest) %>%
  arrange(match(gene_name, txome_genes_interest)) %>%
  .$gene_id %>%
  set_names(txome_genes_interest)

txome_quant_sel <- txome_genes_sel %>%
  imap(\(gene_id, gene_name) {
    txome_quant_small %>%
      plotCounts(
        gene = gene_id,
        intgroup = c("condition_ID"),
        returnData = TRUE
      ) %>%
      droplevels() %>%
      rownames_to_column(var = "name") %>%
      mutate(gene = gene_name)
  }) %>%
  rbindlist() %>%
  as.data.frame() %>%
  mutate(
    gene = dplyr::case_match(
      gene,
      "POU5F1" ~ "OCT3/4",
      .default = gene
    )
  )

# Merge into one table
merged_quant_sel <- rkodata_quant_sel %>%
  group_by(gene) %>%
  left_join(
    mapping,
    by = join_by(name == names_rkodata)
  ) %>%
  left_join(
    txome_quant_sel,
    by = join_by(
      names_txome == name,
      gene == gene
    )
  )

# Create correlation plot
rkodata_plots_corr <- rkodata_genes_interest %>%
  map(\(gene_sel) {
    # Filter data
    data <- merged_quant_sel %>%
      filter(gene == gene_sel)

    # Plot
    plot <- data %>%
      ggplot(
        aes(
          x = dct_neg_mean,
          y = log(count)
        )
      ) +
      geom_point(
        aes(colour = condition_id)
      ) +
      geom_smooth(
        aes(colour = "black"),
        method = "lm",
        se = FALSE
      ) +
      stat_cor(
        aes(label = ..r.label..),
        r.accuracy = 0.01
      ) +
      scale_color_manual(
        values = c(
          palette_merge[c(1, 2, 3, 4, 5, 6)],
          "#000000"
        )
      ) +
      theme_classic() +
      theme(legend.position = "none") +
      ggtitle(
        label = gene_sel %>%
          str_replace_all("_", " ") %>%
          str_wrap(
            width = 20,
            whitespace_only = TRUE
          )
      ) +
      xlab("Norm. log expression") +
      ylab("Norm. counts")

    # Return data
    return(plot)
  }) %>%
  set_names(rkodata_genes_interest)

saveRDS(
  rkodata_plots_corr,
  here::here(
    "output",
    "plots_box",
    "data",
    "rkodata_plots_corr.rds"
  )
)
