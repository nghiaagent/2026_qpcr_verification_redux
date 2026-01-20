# Import packages
library("EnhancedVolcano")
library("limma")
library("tidyverse")

# Load data
fit_contrasts <- readRDS(
  file = here(
    "output",
    "fit_contrasts.RDS"
  )
)

# Set plot specs
cutoff_logfc <- 4
cutoff_padj <- 1e-15
alpha <- 0.05

# Create list of topTables for boxplots and volcano plots
contrasts <- seq_len(ncol(fit_contrasts$contrasts)) %>%
  setNames(colnames(fit_contrasts$contrasts))

## Unclipped for boxplots
list_toptable <- contrasts %>%
  map(\(coef) {
    fit_contrasts %>%
      topTable(
        coef = coef,
        n = Inf,
        sort.by = "none"
      )
  })

## Clipped for volcano plots
list_toptable_clip <- contrasts %>%
  map(\(coef) {
    fit_contrasts %>%
      topTable(
        coef = coef,
        n = Inf,
        sort.by = "none"
      ) %>%
      clip_results(
        cutoff_logfc = cutoff_logfc,
        cutoff_padj = cutoff_padj,
        alpha = alpha
      )
  })

# Draw volcano plots
plots_volcano <- list_toptable_clip %>%
  imap(\(toptable, title) {
    EnhancedVolcano(
      toptable,
      lab = toptable$GENENAME,
      x = "logFC",
      y = "adj.P.Val",
      xlim = c(cutoff_logfc * -1, cutoff_logfc),
      ylim = c(0, -log10(cutoff_padj)),
      ylab = bquote(~ -Log[10] ~ "adjusted p-value"),
      axisLabSize = 8,
      title = title,
      titleLabSize = 8,
      subtitle = NULL,
      caption = NULL,
      pCutoff = 0.05,
      FCcutoff = 1,
      pointSize = toptable$volcano_size,
      labSize = 2,
      boxedLabels = FALSE,
      shapeCustom = toptable$volcano_shape,
      legendPosition = "none",
      drawConnectors = TRUE,
      widthConnectors = 0.1,
      colConnectors = "grey60",
      arrowheads = FALSE,
      gridlines.major = FALSE,
      gridlines.minor = FALSE
    )
  })
