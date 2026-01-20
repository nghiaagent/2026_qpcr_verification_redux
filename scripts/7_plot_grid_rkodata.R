# Import packages
library("DESeq2")
library("tidyverse")

# Load data
rkodata_plots_goi <- readRDS(
  here::here(
    "output",
    "plots_box",
    "data",
    "rkodata_plots_goi.rds"
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

rkodata_plots_corr <- readRDS(
  here::here(
    "output",
    "plots_box",
    "data",
    "rkodata_plots_corr.rds"
  )
)

## Get legend
rkodata_verif_legend <- get_legend(
  rkodata_plots_goi[[1]] +
    guides(color = guide_legend(ncol = 6)) +
    theme(legend.title = element_blank())
)

# Process plots
## Remove legends as appropriate
rkodata_plots_goi <- rkodata_plots_goi %>%
  map(\(plot) {
    plot +
      theme(
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  })

txome_plots_goi <- txome_plots_goi %>%
  map(\(plot) {
    plot +
      theme(
        plot.title = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  })

rkodata_plots_corr <- rkodata_plots_corr %>%
  map(\(plot) {
    plot +
      theme(
        plot.title = element_blank(),
        legend.position = "none"
      )
  })

# Build grid of boxplots
## Build grids
## HSPG core proteins: Indices 1 - 12
## HSPG enzymes: Indices 13 - 22
## MSC markers: Indices 23 - 34
rkodata_merged_plots <- list(
  rkodata_plots_goi,
  txome_plots_goi,
  rkodata_plots_corr
) %>%
  pmap(\(plot_1, plot_2, plot_3) {
    plot <- plot_1 | plot_2 | plot_3

    return(plot)
  })

rkodata_verif_plots <- list(
  hspg_core = c(1:6),
  hspg_enzymes = c(7:11),
  msc_markers = c(12:17)
) %>%
  map(\(index) {
    rkodata_merged_plots[index] %>%
      wrap_plots(
        ncol = 2
      ) /
      rkodata_verif_legend +
      plot_layout(heights = c(20, 1))
  })

# Save data
imap(
  rkodata_verif_plots,
  \(plot, name) {
    ggsave(
      filename = str_c("rkodata_", name, ".png"),
      plot = plot,
      path = here::here(
        "output",
        "plots_box"
      ),
      scale = 0.95,
      width = 13,
      height = 8,
      units = "in",
      dpi = 144
    )
  }
)
