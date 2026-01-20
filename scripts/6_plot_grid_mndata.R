# Declare location
here::i_am("scripts/6_plot_grid_mndata.R")

# Import packages
library("cowplot")
library("patchwork")
library("tidyverse")
library("vctrs")

# Load data
mndata_plots <- readRDS(
  here::here(
    "output",
    "plots_box",
    "data",
    "mndata_plots.rds"
  )
)

txome_plots_goi <- readRDS(
  here::here(
    "output",
    "plots_box",
    "data",
    "txome_plots_goi.rds"
  )
)

mndata_plots_corr <- readRDS(
  here::here(
    "output",
    "plots_box",
    "data",
    "mndata_plots_corr.rds"
  )
)

# Plot as grid
## txome vs MN qpcr
### Collect relevant legend
mndata_plots_hmsc_legend <- get_legend(
  mndata_plots[["hMSC_response"]][[1]] +
    guides(
      colour = guide_legend(ncol = 4),
      fill = "none",
      shape = "none"
    ) +
    theme(legend.title = element_blank())
)

### Create plot grid
mndata_plots_hmsc <- c(
  mndata_plots[["hMSC_response"]] %>%
    map(\(plot) {
      plot +
        theme(
          legend.position = "none",
          axis.title.x = element_blank()
        )
    }),
  txome_plots_goi %>%
    map(\(plot) {
      plot +
        theme(
          plot.title = element_blank(),
          legend.position = "none",
          axis.title.x = element_blank()
        )
    }),
  mndata_plots_corr %>%
    map(\(plot) {
      plot +
        theme(
          plot.title = element_blank(),
          legend.position = "none"
        )
    })
) %>%
  wrap_plots(
    ncol = 3,
    nrow = 4,
    byrow = FALSE
  ) /
  mndata_plots_hmsc_legend +
  plot_layout(heights = c(20, 1))

mndata_plots_hmsc_paper <- vec_interleave(
  mndata_plots[["hMSC_response"]] %>%
    map(\(plot) {
      plot +
        theme(
          legend.position = "none",
          axis.title.x = element_blank()
        )
    }),
  txome_plots_goi %>%
    map(\(plot) {
      plot +
        theme(
          plot.title = element_blank(),
          legend.position = "none",
          axis.title.x = element_blank()
        )
    })
) %>%
  wrap_plots(
    ncol = 4,
    nrow = 2,
    byrow = TRUE
  ) /
  mndata_plots_hmsc_legend +
  plot_layout(heights = c(20, 1))

## Cell populations (baseline, response to hep, log2fc)
### Collect legend
mndata_plots_cellpops_legend <- get_legend(
  mndata_plots[["cellpops_response"]][[1]] +
    guides(
      colour = guide_legend(ncol = 6),
      fill = "none",
      shape = "none"
    ) +
    theme(legend.title = element_blank())
)

### Create plot grid
mndata_plots_cellpops <- list(
  mndata_plots[["cellpops_baseline"]] %>%
    map(\(plot) {
      plot +
        theme(
          legend.position = "none",
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank()
        )
    }),
  mndata_plots[["cellpops_response"]] %>%
    map(\(plot) {
      plot +
        theme(
          plot.title = element_blank(),
          legend.position = "none",
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank()
        )
    }),
  mndata_plots[["cellpops_response_lfc"]] %>%
    map(\(plot) {
      plot +
        theme(
          plot.title = element_blank(),
          legend.position = "none",
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank()
        )
    })
) %>%
  # Group baseline, response to hep, log2fc plots of each gene
  pmap(\(plot_1, plot_2, plot_3) {
    plot <- plot_1 | plot_2 | plot_3

    return(plot)
  }) %>%
  # Create grid
  wrap_plots(ncol = 2) %>%
  # Add legend
  wrap_plots(
    mndata_plots_cellpops_legend,
    ncol = 1,
    heights = c(20, 1)
  )

mndata_plot_cellpops_paper <- mndata_plots[["cellpops_response"]] %>%
  map(\(plot) {
    plot +
      theme(
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()
      )
  }) %>%
  # Create grid
  wrap_plots(ncol = 2) %>%
  # Add legend
  wrap_plots(
    mndata_plots_cellpops_legend,
    ncol = 1,
    heights = c(20, 2)
  )

mndata_plots_cellpops_supplementary <- list(
  mndata_plots[["cellpops_baseline"]] %>%
    map(\(plot) {
      plot +
        theme(
          legend.position = "none",
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank()
        )
    }),
  mndata_plots[["cellpops_response_lfc"]] %>%
    map(\(plot) {
      plot +
        theme(
          plot.title = element_blank(),
          legend.position = "none",
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank()
        )
    })
) %>%
  # Group baseline, response to hep, log2fc plots of each gene
  pmap(\(plot_1, plot_2) {
    plot <- plot_1 | plot_2

    return(plot)
  }) %>%
  # Create grid
  wrap_plots(ncol = 2) %>%
  # Add legend
  wrap_plots(
    mndata_plots_cellpops_legend,
    ncol = 1,
    heights = c(20, 1)
  )

## Cell populations response (relexp)
### Collect legend
mndata_plots_cellpops_response_linexp_legend <- get_legend(
  mndata_plots[["cellpops_response_linexp"]][[1]] +
    guides(fill = guide_legend(ncol = 6, nrow = 2)) +
    theme(legend.title = element_blank())
)

### Collect plot grid
mndata_plots_cellpops_response_linexp <- mndata_plots[[
  "cellpops_response_linexp"
]] %>%
  map(\(plot) {
    plot +
      theme(
        legend.position = "none"
      )
  }) %>%
  wrap_plots(
    ncol = 2,
    byrow = FALSE
  ) /
  mndata_plots_cellpops_response_linexp_legend +
  plot_layout(heights = c(20, 1.5))

# Export plots
## txome vs MN qpcr
ggsave(
  filename = "mndata_hmsc.png",
  plot = mndata_plots_hmsc,
  path = here::here(
    "output",
    "plots_box"
  ),
  scale = 0.8,
  width = 8,
  height = 11,
  units = "in",
  dpi = 144
)

ggsave(
  filename = "mndata_hmsc_paper.png",
  plot = mndata_plots_hmsc_paper,
  path = here::here(
    "output",
    "plots_box"
  ),
  scale = 0.8,
  width = 12,
  height = 6,
  units = "in",
  dpi = 144
)

## Cell populations baseline
ggsave(
  filename = "mndata_cellpops.png",
  plot = mndata_plots_cellpops,
  path = here::here(
    "output",
    "plots_box"
  ),
  scale = 0.9,
  width = 14,
  height = 8,
  units = "in",
  dpi = 144
)

ggsave(
  filename = "mndata_cellpops_paper.png",
  plot = mndata_plot_cellpops_paper,
  path = here::here(
    "output",
    "plots_box"
  ),
  scale = 0.6,
  width = 12,
  height = 12,
  units = "in",
  dpi = 144
)

ggsave(
  filename = "mndata_cellpops_supplementary.png",
  plot = mndata_plots_cellpops_supplementary,
  path = here::here(
    "output",
    "plots_box"
  ),
  scale = 0.8,
  width = 12,
  height = 10,
  units = "in",
  dpi = 144
)

## Cell populations response (linear)
ggsave(
  filename = "mndata_cellpops_response_linexp.png",
  plot = mndata_plots_cellpops_response_linexp,
  path = here::here(
    "output",
    "plots_box"
  ),
  scale = 0.8,
  width = 8,
  height = 8,
  units = "in",
  dpi = 144
)
