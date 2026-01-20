# Import packages
library("DESeq2")
library("magrittr")
library("tidyverse")

# Define conditions of interest
txome_conditions_interest <- c(
  "P5D5Untreated",
  "P13D5Untreated",
  "P5D5Treated",
  "P13D5Treated"
)
# Load data
txome_quant_deseq2_batchcor <- readRDS(
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

# Subset data based on conditions and genes
txome_quant_small <- txome_quant_deseq2_batchcor %$%
  .[, colData(.)$condition_ID %in% txome_conditions_interest] %$%
  .[, colData(.)$cell_line %in% c("hMSC-21558")]

txome_genes_sel <- txome_quant_small %>%
  rowRanges() %>%
  as.data.frame() %>%
  filter(symbol %in% mndata_genes) %>%
  arrange(match(gene_name, mndata_genes))

## Readjust DESeq results to include only GOIs
txome_results_subset <- txome_results %>%
  map(\(x) {
    x <- x[rownames(x) %in% txome_genes_sel$gene_id, ]
    x$padj <- p.adjust(x$pvalue, method = "BH")
    return(x)
  })

# Define functions for drawing boxplots
## Draw boxplot hMSC response to heparin (txome)
txome_plot_goi <- function(gene_sel, gene_name) {
  # Get table containing normalised gene counts
  gene_counts <- txome_quant_small %>%
    plotCounts(
      gene = gene_sel,
      intgroup = c("condition_ID", "Treatment"),
      returnData = TRUE
    ) %>%
    droplevels()

  # Get position to start drawing signif bars
  y_position <- log10(max(gene_counts$count) * 1.1)

  # Define position of bars and numbers manually

  # Plot data
  plot <- gene_counts %>%
    ggplot(aes(
      x = condition_ID,
      y = count
    )) +
    geom_boxplot(
      aes(
        colour = condition_ID
      ),
      outliers = FALSE
    ) +
    geom_jitter(
      aes(
        fill = condition_ID,
        shape = Treatment
      ),
      colour = "black",
      size = 1.2,
      stroke = 0.5,
      alpha = 0.5
    ) +
    geom_signif(
      comparisons = list(
        c("P5D5Untreated", "P5D5Treated"),
        c("P5D5Untreated", "P13D5Untreated"),
        c("P13D5Untreated", "P13D5Treated")
      ),
      annotation = c(2, 21, 6) %>%
        map(\(coef) txome_results_subset[[coef]][gene_sel, ]$padj) %>%
        unlist() %>%
        stars_pval(),
      y_position = y_position %>%
        rep(3),
      textsize = 3,
      step_increase = 0.2
    ) +
    scale_colour_manual(
      values = palette_merge[c(1, 2, 5, 6)],
      labels = c(
        "P5D5Untreated" = "P+5 D5 Control",
        "P5D5Treated" = "P+5 D5 Heparin",
        "P13D5Untreated" = "P+13 D5 Control",
        "P13D5Treated" = "P+13 D5 Heparin"
      )
    ) +
    scale_fill_manual(
      values = palette_merge[c(1, 2, 5, 6)],
      labels = c(
        "P5D5Untreated" = "P+5 D5 Control",
        "P5D5Treated" = "P+5 D5 Heparin",
        "P13D5Untreated" = "P+13 D5 Control",
        "P13D5Treated" = "P+13 D5 Heparin"
      )
    ) +
    scale_shape_manual(
      values = c(
        "Untreated" = 21,
        "Treated" = 24
      )
    ) +
    scale_x_discrete(
      labels = c(
        "P5D5Untreated" = "P+5 D5 Control",
        "P5D5Treated" = "P+5 D5 Heparin",
        "P13D5Untreated" = "P+13 D5 Control",
        "P13D5Treated" = "P+13 D5 Heparin"
      )
    ) +
    scale_y_log10(expand = expansion(0, 0.15)) +
    theme_classic() +
    ggtitle(label = gene_name) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    ylab("Norm. counts")

  # Return object
  return(plot)
}

# Draw plots
## Boxplots of norm counts
txome_plots_goi <- txome_genes_sel %$%
  map2(
    .$gene_id,
    .$gene_name,
    \(x, y) txome_plot_goi(x, y)
  ) %>%
  set_names(txome_genes_sel$gene_name)

# Save data
saveRDS(
  txome_plots_goi,
  here::here(
    "output",
    "plots_box",
    "data",
    "txome_plots_goi.rds"
  )
)
