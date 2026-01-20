# Import packages
library("data.table")
library("DESeq2")
library("tidyverse")

# Define variables
## List of conditions
mndata_conditions_interest <- c(
  "hMSC-21558 P5 D5 Control",
  "hMSC-21558 P5 D5 Heparin",
  "hMSC-21558 P13 D5 Control",
  "hMSC-21558 P13 D5 Heparin"
)
txome_conditions_interest <- c(
  "P5D5Treated",
  "P5D5Untreated",
  "P13D5Untreated",
  "P13D5Treated"
)

## Condition mappings
conditions_interest_map <- c(
  "hMSC-21558 P5 D5 Control" = "P5D5Untreated",
  "hMSC-21558 P5 D5 Heparin" = "P5D5Treated",
  "hMSC-21558 P13 D5 Control" = "P13D5Untreated",
  "hMSC-21558 P13 D5 Heparin" = "P13D5Treated"
)

# Load data
## MN qPCR data
mndata_quant_big <- readRDS(
  file = here(
    "output",
    "quant_big.RDS"
  )
)

## Txome data
txome_quant_deseq2_batchcor <- readRDS(
  file = here::here(
    "input",
    "txome",
    "quant_deseq2_batchcor_nofilter.RDS"
  )
)

# Clean data
# Get table of per biological replicate gene expression

## qPCR data
mndata_quant_sel <- mndata_quant_big %>%
  # Subset data based on genes and condition
  dplyr::filter(sample_id %in% mndata_conditions_interest) %>%
  # Summarise to biol replicate level
  dplyr::group_by(name, gene) %>%
  dplyr::summarise(
    sample_id = sample_id %>%
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
txome_quant_small <- txome_quant_deseq2_batchcor %$%
  .[, colData(.)$condition_ID %in% txome_conditions_interest] %$%
  .[, colData(.)$cell_line %in% c("hMSC-21558")]

### Get gene ids
txome_genes_sel <- txome_quant_small %>%
  rowRanges() %>%
  as.data.frame() %>%
  filter(symbol %in% mndata_genes) %>%
  arrange(match(gene_name, mndata_genes)) %>%
  .$gene_id %>%
  set_names(mndata_genes)

### Get table of counts
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
  rbindlist()

# Merge tables
merged_quant_sel <- full_join(
  mndata_quant_sel,
  txome_quant_sel,
  by = join_by(
    name == name,
    gene == gene
  ),
  suffix = c("qpcr", "txome")
)

# Create correlation plot
mndata_plots_corr <- mndata_genes %>%
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
        aes(colour = condition_ID)
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
          palette_merge[c(1, 2, 5, 6)],
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
  })

# Save data
saveRDS(
  mndata_plots_corr,
  here::here(
    "output",
    "plots_box",
    "data",
    "mndata_plots_corr.rds"
  )
)
