# Import packages
library("DESeq2")
library("tidyverse")

# Plot logFC of txome and Rachel's qPCR data as XY plots

# Set genes for exclusion and GOIs
## RKO data
rkodata_genes_interest <- c(
  "EXT1",
  "EXT2",
  "HS2ST1",
  "NDST1",
  "NDST2",
  "GPC1",
  "GPC4",
  "SDC1",
  "SDC2",
  "SDC3",
  "SDC4",
  "DCN",
  "VCAN",
  "ACTA2",
  "CD44",
  "COL1A1",
  "CTNNB1",
  "FN1",
  "VIM",
  "ENO2",
  "NANOG",
  "OCT3/4",
  "TUBB3"
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
  "EXT1",
  "EXT2",
  "HS2ST1",
  "NDST1",
  "NDST2",
  "GPC1",
  "GPC4",
  "SDC1",
  "SDC2",
  "SDC3",
  "SDC4",
  "DCN",
  "VCAN",
  "ACTA2",
  "CD44",
  "COL1A1",
  "CTNNB1",
  "FN1",
  "VIM",
  "ENO2",
  "NANOG",
  "POU5F1",
  "TUBB3"
)

# Set contrasts of interest
rkodata_contrasts <- c(
  "treatment_P5_D3" = 1,
  "treatment_P7_D3" = 2,
  "treatment_P13_D3" = 3,
  "P7vsP5_D3" = 7,
  "P13vsP5_D3" = 8,
  "P13vsP7_D3" = 9
)
txome_contrasts <- c(
  "Trt_P5_D3",
  "Trt_P7_D3",
  "Trt_P13_D3",
  "P7vsP5_UT_D3",
  "P13vsP5_UT_D3",
  "P13vsP7_UT_D3"
)

# Load data
## RKO data
### Limma fit of heparin data
rkodata_fit_contrasts <- readRDS(
  here::here(
    "input",
    "rkodata",
    "4_fit_contrasts.rds"
  )
) %>%
  # Select heparin contrasts from list
  .[["heparin"]] %>%
  # Filter for GOIs
  .[rkodata_genes_interest, ]

rkodata_results <- map(
  rkodata_contrasts,
  \(coef) {
    topTable(
      rkodata_fit_contrasts,
      coef = coef,
      number = Inf,
      sort.by = "none"
    )
  }
)

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
### Subset DESeqResults based on GOIs
#### Get list of gene ids
txome_genes_interest_id <- txome_quant_deseq2 %>%
  rowRanges() %>%
  as.data.frame() %>%
  filter(symbol %in% txome_genes_interest) %>%
  arrange(match(gene_name, txome_genes_interest))

#### Readjust DESeq results to include only GOIs
txome_results_subset <- txome_results %>%
  .[txome_contrasts] %>%
  map(\(x) {
    x <- x[rownames(x) %in% txome_genes_interest_id$gene_id, ]
    x$padj <- p.adjust(x$pvalue, method = "BH")
    return(x)
  })

# Merge results tables
merged_results <- map2(
  txome_results_subset,
  rkodata_results,
  \(txome, rkodata) {
    cbind(txome, rkodata) %>%
      as.data.frame() %>%
      # Add outcome column
      mutate(
        outcome_txome = case_when(
          padj < 0.05 & log2FoldChange > 0 ~ 1,
          padj < 0.05 & log2FoldChange < 0 ~ -1,
          .default = 0
        ),
        outcome_rkodata = case_when(
          adj.P.Val < 0.05 & logFC > 0 ~ 1,
          adj.P.Val < 0.05 & logFC < 0 ~ -1,
          .default = 0
        )
      )
  }
)

# Draw EnhancedVolcano
plots <- merged_results %>%
  imap(\(toptable, name) {
    toptable %>%
      ggplot(
        aes(
          x = log2FoldChange,
          y = logFC,
          fill = outcome2factor(outcome_txome),
          colour = outcome2factor(outcome_rkodata)
        )
      ) +
      geom_point(
        shape = "circle filled",
        stroke = 1
      ) +
      geom_quadrant_lines(linetype = "dotted") +
      scale_fill_outcome() +
      scale_colour_outcome() +
      scale_x_continuous(limits = symmetric_limits) +
      scale_y_continuous(limits = symmetric_limits) +
      coord_obs_pred() +
      theme_classic() +
      ggtitle(name)
  })

plots_legend <- get_legend(
  plots[[6]] +
    guides(
      colour = guide_legend(ncol = 3),
      fill = guide_legend(ncol = 3)
    ) +
    theme(legend.title = element_blank())
)

plots_out <- plots %>%
  map(\(plot) {
    plot +
      xlab("logFC (txome)") +
      ylab("logFC (RT-qPCR)") +
      theme(legend.position = "none")
  }) %>%
  wrap_plots(
    ncol = 3
  ) /
  plots_legend +
  plot_layout(heights = c(20, 1))

# Save plots
ggsave(
  filename = "rkodata_quadrant.png",
  plot = plots_out,
  path = here::here(
    "output",
    "plots_box"
  ),
  scale = 0.75,
  width = 9,
  height = 7,
  units = "in",
  dpi = 144
)
