### Load data

source("./scripts/analyse_limma.R")

p_load(ggbeeswarm,
       ggpattern,
       ggh4x,
       palettetown)

# Make gene list

genes <- tibble(name = rownames(dge$E))
genes$group <-
  map(genes$name, \(y) names(which(sapply(gene_groups, function(x)
    y %in% x)))) %>%
  factor(levels = names(gene_groups))

genes <- genes %>%
  arrange(group, name)

facet_passage <- c("P5",
                   "P7",
                   "P13")

names(facet_passage) <- c("5",
                          "7",
                          "13")

facet_day <- c("D3",
               "D5")

names(facet_day) <- c("3",
                      "5")

# Set up function

draw_trellis <- function(gene_sel) {
  df_plot <- filter(df_big, Detector == gene_sel)
  
  df_summary <- df_plot %>%
    group_by(passage, day, treatment_group) %>%
    summarise(
      mean = mean(dCt_nexp_log2_offset),
      se = sd(dCt_nexp_log2_offset) / sqrt(length(dCt_nexp_log2_offset))
    )
  
  ggplot() +
    geom_boxplot_pattern(
      data = df_plot,
      aes(
        x = treatment_group,
        y = dCt_nexp_log2_offset,
        pattern_density = treatment_group,
      ),
      outliers = FALSE,
      fill = "white",
      colour = 'black',
      pattern = 'stripe',
      pattern_alpha = 0.4
    ) +
    scale_pattern_density_manual(
      values = c('Ctrl' = 0, '50mM' = 0.7),
      labels = c('Ctrl' = "Control",
                 '50mM' = 'Chlorate')
    ) +
    # geom_errorbar(data = df_summary,
    #               aes(
    #                 x = treatment_group,
    #                 ymin = (mean),
    #                 ymax = (mean + se),
    #               ),
    #               width = 0.4) +
    geom_quasirandom(
      data = df_plot,
      aes(x = treatment_group,
          y = dCt_nexp_log2_offset,
          color = cell_line),
      alpha = 0.8,
      size = 1
    ) +
    labs(
      x = "Treatment",
      y = expression(paste("Relative log-expression (+50)")),
      colour = "Cell population",
      pattern_density = "Treatment"
    ) +
    scale_x_discrete(labels = c("Ctrl" = "Control",
                                "50mM" = "Chlorate")) +
    scale_colour_manual(values = c(palettetown::pokepal(6)[3],
                                   palettetown::pokepal(6)[4])) +
    facet_nested(
      ~ passage + day,
      switch = 'x',
      labeller = labeller(passage = facet_passage,
                          day = facet_day)
    ) +
    theme_minimal(base_size = 22) +
    theme(
      legend.position = "none",
      legend.key.width = unit(0.5, 'in'),
      panel.grid.major.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(
        family = "Georgia",
        size = rel(1.5),
        vjust = 0.3,
        hjust = 1
      ),
      strip.text = element_text(size = rel(1.2))
    )
}

plots <-
  list(
    "GPC1" = draw_trellis('GPC1'),
    "GPC4" = draw_trellis('GPC4'),
    "NDST1" = draw_trellis('NDST1'),
    "NDST2" = draw_trellis('NDST2'),
    "EXT1" = draw_trellis('EXT1'),
    "EXT2" = draw_trellis('EXT2'),
    "GLCE" = draw_trellis('GLCE'),
    "HS2ST1" = draw_trellis('HS2ST1'),
    "HS6ST1" = draw_trellis('HS6ST1'),
    "SDC1" = draw_trellis('SDC1'),
    "SDC2" = draw_trellis('SDC2'),
    "SDC3" = draw_trellis('SDC3'),
    "SDC4" = draw_trellis('SDC4'),
    "NES" = draw_trellis('NES'),
    "TUBB3" = draw_trellis('TUBB3'),
    "NANOG" = draw_trellis('NANOG'),
    "OCT3/4" = draw_trellis('OCT3/4')
  )

map2(plots,
     names(plots),
     \ (x, y)
     ggsave(
       x,
       filename = file.path(".",
                            "plots",
                            "202406_barplots",
                            str_c(y,
                                  ".png",
                                  sep = "")),
       width = 10,
       height = 6,
       scale = 1
     ))