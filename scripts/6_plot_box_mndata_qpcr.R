# Import packages
library("tidyverse")
library("ggpattern")
library("ggsignif")

# Load data
mndata_quant_big <- readRDS(
  file = here(
    "output",
    "quant_big.RDS"
  )
)
mndata_fit_contrasts <- readRDS(
  file = here(
    "output",
    "fit_contrasts.RDS"
  )
)

# Set contrasts of interest
mndata_contrasts_interest_response <- c(
  "Hep_hMSC_P5",
  "Hep_hMSC_P13",
  "Hep_RCX",
  "Hep_RVM",
  "Hep_1321N1",
  "Hep_SH"
)

# Create list of topTables for boxplots and volcano plots
# Unclipped for boxplots
mndata_contrasts <- seq_len(ncol(mndata_fit_contrasts$contrasts)) %>%
  setNames(colnames(mndata_fit_contrasts$contrasts))

mndata_toptable <- mndata_contrasts %>%
  map(\(coef) {
    mndata_fit_contrasts %>%
      topTable(
        coef = coef,
        n = Inf,
        sort.by = "none",
        confint = TRUE
      )
  })

# Define functions for drawing boxplots
## Draw boxplot hMSC response to heparin
mndata_plot_goi_hMSC_response <- function(gene_sel) {
  # Get table containing dCt for hMSC
  quant_small <- mndata_quant_big %>%
    dplyr::filter(
      gene %in% gene_sel,
      cell_line == "hMSC-21558"
    ) %>%
    droplevels()

  # Get position to start drawing signif bars
  y_position <- max(quant_small$dct_nexp_log2) * 0.95

  # Plot data
  plot <- quant_small %>%
    # Define plot x y axes
    ggplot(aes(
      x = sample_id,
      y = dct_nexp_log2
    )) +
    geom_boxplot(
      aes(
        colour = sample_id
      ),
      outliers = FALSE
    ) +
    geom_jitter(
      aes(
        fill = sample_id,
        shape = treatment
      ),
      color = "black",
      size = 1.2,
      stroke = 0.5,
      alpha = 0.5
    ) +
    geom_signif(
      comparisons = list(
        c("hMSC-21558 P5 D5 Control", "hMSC-21558 P5 D5 Heparin"),
        c("hMSC-21558 P5 D5 Control", "hMSC-21558 P13 D5 Control"),
        c("hMSC-21558 P13 D5 Control", "hMSC-21558 P13 D5 Heparin")
      ),
      annotation = c(3, 1, 2) %>%
        map(\(x) mndata_toptable[[x]][gene_sel, "adj.P.Val"]) %>%
        unlist() %>%
        stars_pval(),
      y_position = y_position %>%
        rep(3),
      textsize = 3,
      step_increase = 0.2
    ) +
    scale_y_continuous(expand = expansion(0, 0.4)) +
    theme_classic() +
    ggtitle(
      label = gene_sel %>%
        str_replace_all("_", " ") %>%
        str_wrap(
          width = 20,
          whitespace_only = TRUE
        )
    ) +
    scale_x_discrete(
      labels = c(
        "hMSC-21558 P5 D5 Control" = "P+5 D5 Control",
        "hMSC-21558 P5 D5 Heparin" = "P+5 D5 Heparin",
        "hMSC-21558 P13 D5 Control" = "P+13 D5 Control",
        "hMSC-21558 P13 D5 Heparin" = "P+13 D5 Heparin"
      )
    ) +
    scale_colour_manual(
      values = palette_merge[c(1, 2, 5, 6)],
      labels = c(
        "hMSC-21558 P5 D5 Control" = "P+5 D5 Control",
        "hMSC-21558 P5 D5 Heparin" = "P+5 D5 Heparin",
        "hMSC-21558 P13 D5 Control" = "P+13 D5 Control",
        "hMSC-21558 P13 D5 Heparin" = "P+13 D5 Heparin"
      )
    ) +
    scale_fill_manual(
      values = palette_merge[c(1, 2, 5, 6)],
      labels = c(
        "hMSC-21558 P5 D5 Control" = "P+5 D5 Control",
        "hMSC-21558 P5 D5 Heparin" = "P+5 D5 Heparin",
        "hMSC-21558 P13 D5 Control" = "P+13 D5 Control",
        "hMSC-21558 P13 D5 Heparin" = "P+13 D5 Heparin"
      )
    ) +
    scale_shape_manual(
      values = c(
        "Control" = 21,
        "Heparin" = 24
      )
    ) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    ylab("Norm. log expression")

  # Return object
  return(plot)
}

## Draw boxplot cell type differences in baseline
mndata_plot_cellpops_baseline <- function(gene_sel) {
  # Get table containing dCt for control only
  quant_small <- mndata_quant_big %>%
    dplyr::filter(
      gene %in% gene_sel,
      treatment == "Control"
    ) %>%
    droplevels()

  # Get position to start drawing signif bars
  y_position <- max(quant_small$dct_nexp_log2) * 0.95

  # Plot data
  plot <- quant_small %>%
    # Define plot x y axes
    ggplot(aes(
      x = sample_id,
      y = dct_nexp_log2
    )) +
    # Define color of boxplot and jitter
    geom_boxplot(
      aes(
        colour = sample_id
      ),
      outliers = FALSE
    ) +
    geom_jitter(
      aes(
        fill = sample_id,
        shape = treatment
      ),
      colour = "black",
      size = 1.2,
      stroke = 0.5,
      alpha = 0.5
    ) +
    # Define significance signs and bars
    geom_signif(
      # Bars connecting boxplots
      comparisons = list(
        c("hMSC-21558 P5 D5 Control", "ReNcell CX D5 Control"),
        c("hMSC-21558 P5 D5 Control", "ReNcell VM D5 Control"),
        c("hMSC-21558 P5 D5 Control", "1321N1 D5 Control"),
        c("hMSC-21558 P5 D5 Control", "SH-SY5Y D3 Control"),
        c("ReNcell CX D5 Control", "ReNcell VM D5 Control"),
        c("ReNcell CX D5 Control", "1321N1 D5 Control"),
        c("ReNcell VM D5 Control", "SH-SY5Y D3 Control"),
        c("1321N1 D5 Control", "SH-SY5Y D3 Control")
      ),
      # Extract p-vals using coefficients of contrast matrix
      annotation = c(
        10,
        11,
        12,
        13,
        14,
        15,
        18,
        19
      ) %>%
        map(\(x) mndata_toptable[[x]][gene_sel, "adj.P.Val"]) %>%
        unlist() %>%
        stars_pval(),
      # Define y-position of bars
      y_position = y_position %>%
        rep(8),
      textsize = 3,
      step_increase = 0.2
    ) +
    scale_y_continuous(expand = expansion(0, 3)) +
    theme_classic() +
    ggtitle(
      label = gene_sel %>%
        str_replace_all("_", " ") %>%
        str_wrap(
          width = 20,
          whitespace_only = TRUE
        )
    ) +
    scale_x_discrete(
      labels = c(
        "hMSC-21558 P5 D5 Control" = "hMSC P+5",
        "hMSC-21558 P13 D5 Control" = "hMSC P+13",
        "ReNcell CX D5 Control" = "ReNcell CX",
        "SH-SY5Y D3 Control" = "SH-SY5Y",
        "ReNcell VM D5 Control" = "ReNcell VM",
        "1321N1 D5 Control" = "1321N1"
      )
    ) +
    scale_colour_manual(
      values = palette_merge[c(1, 5, 7, 9, 11, 13)],
      labels = c(
        "hMSC-21558 P5 D5 Control" = "hMSC P+5",
        "hMSC-21558 P13 D5 Control" = "hMSC P+13",
        "ReNcell CX D5 Control" = "ReNcell CX",
        "SH-SY5Y D3 Control" = "SH-SY5Y",
        "ReNcell VM D5 Control" = "ReNcell VM",
        "1321N1 D5 Control" = "1321N1"
      )
    ) +
    scale_fill_manual(
      values = palette_merge[c(1, 5, 7, 9, 11, 13)],
      labels = c(
        "hMSC-21558 P5 D5 Control" = "hMSC P+5",
        "hMSC-21558 P13 D5 Control" = "hMSC P+13",
        "ReNcell CX D5 Control" = "ReNcell CX",
        "SH-SY5Y D3 Control" = "SH-SY5Y",
        "ReNcell VM D5 Control" = "ReNcell VM",
        "1321N1 D5 Control" = "1321N1"
      )
    ) +
    scale_shape_manual(
      values = c(
        "Control" = 21,
        "Heparin" = 24
      )
    ) +
    ylab("Norm. log expression")

  # Return object
  return(plot)
}

## Draw boxplot cell type differences in response (log exprsn)
mndata_plot_cellpops_response <- function(gene_sel) {
  # Get table containing dCt for control only
  quant_small <- mndata_quant_big %>%
    dplyr::filter(gene %in% gene_sel) %>%
    droplevels()

  # Get position to start drawing signif bars
  y_position <- max(quant_small$dct_nexp_log2) * 0.95

  # Plot data
  plot <- quant_small %>%
    # Define plot x y axes
    ggplot(aes(
      x = sample_id,
      y = dct_nexp_log2
    )) +
    geom_boxplot(
      aes(
        colour = sample_id
      ),
      outliers = FALSE
    ) +
    geom_jitter(
      aes(
        fill = sample_id,
        shape = treatment
      ),
      colour = "black",
      size = 1.2,
      stroke = 0.5,
      alpha = 0.5
    ) +
    geom_signif(
      comparisons = list(
        c("hMSC-21558 P5 D5 Control", "hMSC-21558 P5 D5 Heparin"),
        c("hMSC-21558 P13 D5 Control", "hMSC-21558 P13 D5 Heparin"),
        c("ReNcell CX D5 Control", "ReNcell CX D5 Heparin"),
        c("ReNcell VM D5 Control", "ReNcell VM D5 Heparin"),
        c("1321N1 D5 Control", "1321N1 D5 Heparin"),
        c("SH-SY5Y D3 Control", "SH-SY5Y D3 Heparin")
      ),
      annotation = c(3, 4, 6, 7, 8, 9) %>%
        map(\(x) mndata_toptable[[x]][gene_sel, "adj.P.Val"]) %>%
        unlist() %>%
        stars_pval(),
      y_position = y_position %>%
        rep(6),
      textsize = 3
    ) +
    scale_y_continuous(expand = expansion(0, 3)) +
    theme_classic() +
    ggtitle(
      label = gene_sel %>%
        str_replace_all("_", " ") %>%
        str_wrap(
          width = 20,
          whitespace_only = TRUE
        )
    ) +
    scale_x_discrete(
      labels = c(
        "hMSC-21558 P5 D5 Control" = "hMSC P+5 Control",
        "hMSC-21558 P5 D5 Heparin" = "hMSC P+5 Heparin",
        "hMSC-21558 P13 D5 Control" = "hMSC P+13 Control",
        "hMSC-21558 P13 D5 Heparin" = "hMSC P+13 Heparin",
        "ReNcell CX D5 Control" = "ReNcell CX Control",
        "ReNcell CX D5 Heparin" = "ReNcell CX Heparin",
        "ReNcell VM D5 Control" = "ReNcell VM Control",
        "ReNcell VM D5 Heparin" = "ReNcell VM Heparin",
        "1321N1 D5 Control" = "1321N1 Control",
        "1321N1 D5 Heparin" = "1321N1 Heparin",
        "SH-SY5Y D3 Control" = "SH-SY5Y Control",
        "SH-SY5Y D3 Heparin" = "SH-SY5Y Heparin"
      )
    ) +
    scale_colour_manual(
      values = palette_merge[c(1:2, 5:14)],
      labels = c(
        "hMSC-21558 P5 D5 Control" = "hMSC P+5 Control",
        "hMSC-21558 P5 D5 Heparin" = "hMSC P+5 Heparin",
        "hMSC-21558 P13 D5 Control" = "hMSC P+13 Control",
        "hMSC-21558 P13 D5 Heparin" = "hMSC P+13 Heparin",
        "ReNcell CX D5 Control" = "ReNcell CX Control",
        "ReNcell CX D5 Heparin" = "ReNcell CX Heparin",
        "ReNcell VM D5 Control" = "ReNcell VM Control",
        "ReNcell VM D5 Heparin" = "ReNcell VM Heparin",
        "1321N1 D5 Control" = "1321N1 Control",
        "1321N1 D5 Heparin" = "1321N1 Heparin",
        "SH-SY5Y D3 Control" = "SH-SY5Y Control",
        "SH-SY5Y D3 Heparin" = "SH-SY5Y Heparin"
      )
    ) +
    scale_fill_manual(
      values = palette_merge[c(1:2, 5:14)],
      labels = c(
        "hMSC-21558 P5 D5 Control" = "hMSC P+5 Control",
        "hMSC-21558 P5 D5 Heparin" = "hMSC P+5 Heparin",
        "hMSC-21558 P13 D5 Control" = "hMSC P+13 Control",
        "hMSC-21558 P13 D5 Heparin" = "hMSC P+13 Heparin",
        "ReNcell CX D5 Control" = "ReNcell CX Control",
        "ReNcell CX D5 Heparin" = "ReNcell CX Heparin",
        "ReNcell VM D5 Control" = "ReNcell VM Control",
        "ReNcell VM D5 Heparin" = "ReNcell VM Heparin",
        "1321N1 D5 Control" = "1321N1 Control",
        "1321N1 D5 Heparin" = "1321N1 Heparin",
        "SH-SY5Y D3 Control" = "SH-SY5Y Control",
        "SH-SY5Y D3 Heparin" = "SH-SY5Y Heparin"
      )
    ) +
    scale_shape_manual(
      values = c(
        "Control" = 21,
        "Heparin" = 24
      )
    ) +
    scale_pattern_manual(values = c("none", "stripe")) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    ylab("Norm. log expression")

  # Return object
  return(plot)
}

## Draw boxplot cell type differences in response (logFC)
mndata_plot_cellpops_response_lfc <- function(gene_sel) {
  # Get table containing logFC
  toptable_sel <- mndata_toptable %>%
    imap(\(toptable, name) {
      toptable[gene_sel, ] %>%
        mutate(contrast = name)
    }) %>%
    bind_rows() %>%
    mutate(standard_error = (CI.R - CI.L) / 3.92) %>%
    filter(contrast %in% mndata_contrasts_interest_response) %>%
    mutate(
      contrast = contrast %>%
        factor(levels = mndata_contrasts_interest_response)
    )

  # Get position to start drawing signif bars
  y_position <- max(toptable_sel$logFC + toptable_sel$standard_error)

  ylim <- c(
    max(toptable_sel$logFC + toptable_sel$standard_error) %>%
      abs(),
    min(toptable_sel$logFC - toptable_sel$standard_error) %>%
      abs()
  ) %>%
    max()

  # PLot data
  plot <- toptable_sel %>%
    # Define plot x y axes
    ggplot(aes(x = contrast, y = logFC)) +
    geom_col(
      aes(fill = contrast),
      colour = "grey30"
    ) +
    geom_errorbar(
      aes(
        ymin = logFC - standard_error,
        ymax = logFC + standard_error
      ),
      colour = "grey30",
      width = 0.4,
      alpha = 0.9
    ) +
    geom_signif(
      # Select bars connecting boxplots
      comparisons = list(
        c("Hep_hMSC_P5", "Hep_hMSC_P13"),
        c("Hep_hMSC_P5", "Hep_RCX"),
        c("Hep_hMSC_P5", "Hep_RVM"),
        c("Hep_hMSC_P5", "Hep_1321N1"),
        c("Hep_hMSC_P5", "Hep_SH"),
        c("Hep_RCX", "Hep_RVM"),
        c("Hep_RCX", "Hep_1321N1"),
        c("Hep_RVM", "Hep_SH"),
        c("Hep_1321N1", "Hep_SH")
      ),
      # Extract padj using coefficients
      annotation = c(
        5,
        20,
        21,
        22,
        23,
        24,
        25,
        28,
        29
      ) %>%
        map(\(x) mndata_toptable[[x]][gene_sel, "adj.P.Val"]) %>%
        unlist() %>%
        stars_pval(),
      # Define y-position of bars
      y_position = y_position %>%
        rep(9),
      textsize = 2.4,
      step_increase = 0.25
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed"
    ) +
    scale_y_continuous(
      expand = expansion(0, 3)
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
    scale_x_discrete(
      labels = c(
        "Hep_hMSC_P5" = "hMSC P+5",
        "Hep_hMSC_P13" = "hMSC P+13",
        "Hep_RCX" = "ReNcell CX",
        "Hep_SH" = "SH-SY5Y",
        "Hep_RVM" = "ReNcell VM",
        "Hep_1321N1" = "1321N1"
      )
    ) +
    scale_fill_manual(
      values = palette_merge[c(1, 5, 7, 9, 11, 13)],
      labels = c(
        "Hep_hMSC_P5" = "hMSC P+5",
        "Hep_hMSC_P13" = "hMSC P+13",
        "Hep_RCX" = "ReNcell CX",
        "Hep_SH" = "SH-SY5Y",
        "Hep_RVM" = "ReNcell VM",
        "Hep_1321N1" = "1321N1"
      )
    ) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    ) +
    ylab("log2 fold change")

  # Return object
  return(plot)
}

## Draw boxplot cell type differences in response (linexp)
mndata_plot_cellpops_response_linexp <- function(gene_sel) {
  # Get normalisation factor for dct
  # i.e. Get geomean of hMSC P5 Ctrl, divide all dct by geomean
  factor_norm <- mndata_quant_big %>%
    dplyr::filter(
      gene %in% gene_sel,
      sample_id == "hMSC-21558 P5 D5 Control"
    ) %>%
    .$dct_nexp %>%
    log() %>%
    mean() %>%
    exp()

  # Get table containing dCt for control only
  quant_small <- mndata_quant_big %>%
    dplyr::filter(gene %in% gene_sel) %>%
    droplevels() %>%
    mutate(dct_nexp = dct_nexp / factor_norm)

  # Get table of sumnmary stats
  quant_summary <- quant_small %>%
    dplyr::group_by(sample_id) %>%
    summarise(
      mean = dct_nexp %>%
        log() %>%
        mean() %>%
        exp(),
      se = dct_nexp %>%
        std.error()
    )

  # Get position to start drawing signif bars
  y_position <- max(quant_summary$mean + quant_summary$se)

  # Plot data
  plot <- quant_summary %>%
    # Define plot x y axes
    ggplot(aes(x = sample_id, y = mean)) +
    geom_col(aes(fill = sample_id)) +
    geom_errorbar(
      aes(
        ymin = mean - se,
        ymax = mean + se
      ),
      width = 0.4,
      alpha = 0.9
    ) +
    geom_signif(
      comparisons = list(
        c("hMSC-21558 P5 D5 Control", "hMSC-21558 P5 D5 Heparin"),
        c("hMSC-21558 P13 D5 Control", "hMSC-21558 P13 D5 Heparin"),
        c("ReNcell CX D5 Control", "ReNcell CX D5 Heparin"),
        c("ReNcell VM D5 Control", "ReNcell VM D5 Heparin"),
        c("1321N1 D5 Control", "1321N1 D5 Heparin"),
        c("SH-SY5Y D3 Control", "SH-SY5Y D3 Heparin")
      ),
      annotation = c(3, 4, 6, 7, 8, 9) %>%
        map(\(x) mndata_toptable[[x]][gene_sel, "adj.P.Val"]) %>%
        unlist() %>%
        stars_pval(),
      y_position = y_position %>%
        rep(6),
      textsize = 3
    ) +
    scale_y_continuous(
      expand = expansion(
        mult = c(0, 0.1),
        add = 0
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
    scale_fill_manual(
      values = palette_merge[c(1:2, 5:14)],
      labels = c(
        "hMSC-21558 P5 D5 Control" = "hMSC P+5 Control",
        "hMSC-21558 P5 D5 Heparin" = "hMSC P+5 Heparin",
        "hMSC-21558 P13 D5 Control" = "hMSC P+13 Control",
        "hMSC-21558 P13 D5 Heparin" = "hMSC P+13 Heparin",
        "ReNcell CX D5 Control" = "ReNcell CX Control",
        "ReNcell CX D5 Heparin" = "ReNcell CX Heparin",
        "ReNcell VM D5 Control" = "ReNcell VM Control",
        "ReNcell VM D5 Heparin" = "ReNcell VM Heparin",
        "1321N1 D5 Control" = "1321N1 Control",
        "1321N1 D5 Heparin" = "1321N1 Heparin",
        "SH-SY5Y D3 Control" = "SH-SY5Y Control",
        "SH-SY5Y D3 Heparin" = "SH-SY5Y Heparin"
      )
    ) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    ) +
    ylab("Relative expression")

  # Return object
  return(plot)
}

# Draw list of boxplots
mndata_plots <- list(
  hMSC_response = mndata_genes %>%
    map(\(gene) mndata_plot_goi_hMSC_response(gene)),
  cellpops_baseline = mndata_genes %>%
    map(\(gene) mndata_plot_cellpops_baseline(gene)),
  cellpops_response = mndata_genes %>%
    map(\(gene) mndata_plot_cellpops_response(gene)),
  cellpops_response_lfc = mndata_genes %>%
    map(\(gene) mndata_plot_cellpops_response_lfc(gene)),
  cellpops_response_linexp = mndata_genes %>%
    map(\(gene) mndata_plot_cellpops_response_linexp(gene))
) %>%
  map(\(plotlist) set_names(plotlist, mndata_genes))

# Save data
saveRDS(
  mndata_plots,
  here::here(
    "output",
    "plots_box",
    "data",
    "mndata_plots.rds"
  )
)
