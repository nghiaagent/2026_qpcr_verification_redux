# This script draws up some plots for the CGPH mid year symposium

# Load data

# Get packages
if (!require("pacman"))
  install.packages("pacman")
library(pacman)
pacman::p_load(tidyverse,
               statmod,
               EnhancedVolcano,
               lme4,
               sjPlot,
               lmerTest,
               emmeans,
               afex,
               cowplot)

# Set plot theme

pacman::p_load_gh("Mikata-Project/ggthemr")
ggthemr('pale')
# Load data

df_big <- readRDS("./data/raw_cleaned/df_big.rds")

# List excluded genes

genes_exclude <- c(
  "ALB",
  "COL2A1",
  "FGFR2",
  "FOXA2",
  "GPC2",
  "GPC3",
  "GPC5",
  "HS6ST2",
  "HS6ST3",
  "MAP2",
  "MSI1",
  "NDST3",
  "NDST4",
  "NEUROG2",
  "SOX1"
)

# Data selection
## Look into D3 expression
## Chl treatment
## Modify passages to assist with labelling

df_big <- df_big %>%
  filter(!Detector %in% genes_exclude) %>%
  filter(treatment == "Chl") %>%
  filter(day == "3")

df_big$passage <- factor(df_big$passage,
                         levels = c("5", "7", "13"),
                         labels = c("Phase A", "Phase B", "Phase C"))


# Set up plotting function (by passage)

report_plot <- function(fit) {
  
  by_treatment <- afex_plot(
    fit,
    x = "treatment_group",
    mapping = c("fill", "colour"),
    trace = "cell_line",
    panel = c("passage"),
    error = "model",
    dodge = 0.7,
    line_arg = list(linewidth = 1),
    data_geom = list(ggplot2::geom_boxplot,
                     ggbeeswarm::geom_quasirandom),
    data_arg = list(list(width = 0.5,
                         linetype = 1,
                         alpha = 0.6),
                    list(
                      dodge.width = 0.7,
                      size = 1,
                      alpha = 0.9
                    )),
    id = "ID"
  ) +
    labs(x = NULL,
         y = "Relative expression (Log2)") +
    scale_x_discrete(labels = c('Control', 'Chlorate')) +
    theme(legend.position = "none") +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(by_treatment)
  
}


report_plot_combined <- function(fit) {
  
  by_treatment <- afex_plot(
    fit,
    x = "treatment_group",
    mapping = c("fill", "colour"),
    trace = "treatment_group",
    panel = c("passage"),
    error = "model",
    dodge = 0.7,
    line_arg = list(linewidth = 1),
    data_geom = list(ggplot2::geom_boxplot,
                     ggbeeswarm::geom_quasirandom),
    data_arg = list(list(width = 0.5,
                         linetype = 1,
                         alpha = 0.6),
                    list(
                      dodge.width = 0.7,
                      size = 1,
                      alpha = 0.9
                    )),
    id = "ID"
  ) +
    labs(x = NULL,
         y = "Relative expression (Log2)") +
    scale_x_discrete(labels = c('Control', 'Chlorate')) +
    theme(legend.position = "none") +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(by_treatment)
  
}

# Panel 1: HSPG core proteins

## Plot 1: AGG

df_AGG <- filter(df_big,
                 Detector == "AGG")

# Fit

fit_AGG <-
  mixed(
    dCt_nexp_log2 ~ treatment_group * passage + cell_line + (1 |
                                                               Plate:treatment_group:biological_replicate) + (1 | Plate),
    data = df_AGG
  )

# Plot

plot_AGG <- report_plot(fit_AGG) +
  ggtitle("AGG")

plot_AGG_combined <- report_plot_combined(fit_AGG) +
  ggtitle("AGG")

## Plot 2: GPC1

df_GPC1 <- filter(df_big,
                  Detector == "GPC1")

# Fit

fit_GPC1 <-
  mixed(
    dCt_nexp_log2 ~ treatment_group * passage + cell_line + (1 |
                                                               Plate:treatment_group:biological_replicate) + (1 | Plate),
    data = df_GPC1
  )

# Plot

plot_GPC1 <- report_plot(fit_GPC1) +
  ggtitle("GPC1")

plot_GPC1_combined <- report_plot_combined(fit_GPC1) +
  ggtitle("GPC1")

## Plot 3: GPC4

df_GPC4 <- filter(df_big,
                  Detector == "GPC4")

# Fit

fit_GPC4 <-
  mixed(
    dCt_nexp_log2 ~ treatment_group * passage + cell_line + (1 |
                                                               Plate:treatment_group:biological_replicate) + (1 | Plate),
    data = df_GPC4
  )

# Plot

plot_GPC4 <- report_plot(fit_GPC4) +
  ggtitle("GPC4")
  

plot_GPC4_combined <- report_plot_combined(fit_GPC4) +
  ggtitle("GPC4")



## Plot 4: NDST1

df_NDST1 <- filter(df_big,
                  Detector == "NDST1")

# Fit

fit_NDST1 <-
  mixed(
    dCt_nexp_log2 ~ treatment_group * passage + cell_line + (1 |
                                                               Plate:treatment_group:biological_replicate) + (1 | Plate),
    data = df_NDST1
  )

# Plot

plot_NDST1 <- report_plot(fit_NDST1) +
  ggtitle("NDST1")


plot_NDST1_combined <- report_plot_combined(fit_NDST1) +
  ggtitle("NDST1")


## Plot 5: NDST2

df_NDST2 <- filter(df_big,
                  Detector == "NDST2")

# Fit

fit_NDST2 <-
  mixed(
    dCt_nexp_log2 ~ treatment_group * passage + cell_line + (1 |
                                                               Plate:treatment_group:biological_replicate) + (1 | Plate),
    data = df_NDST2
  )

# Plot

plot_NDST2 <- report_plot(fit_NDST2) +
  ggtitle("NDST2")


plot_NDST2_combined <- report_plot_combined(fit_NDST2) +
  ggtitle("NDST2")


## Plot 6: HS2ST1

df_HS2ST1 <- filter(df_big,
                  Detector == "HS2ST1")

# Fit

fit_HS2ST1 <-
  mixed(
    dCt_nexp_log2 ~ treatment_group * passage + cell_line + (1 |
                                                               Plate:treatment_group:biological_replicate) + (1 | Plate),
    data = df_HS2ST1
  )

# Plot

plot_HS2ST1 <- report_plot(fit_HS2ST1) +
  ggtitle("HS2ST1")


plot_HS2ST1_combined <- report_plot_combined(fit_HS2ST1) +
  ggtitle("HS2ST1")


## Plot 7: HS6ST1

df_HS6ST1 <- filter(df_big,
                  Detector == "HS6ST1")

# Fit

fit_HS6ST1 <-
  mixed(
    dCt_nexp_log2 ~ treatment_group * passage + cell_line + (1 |
                                                               Plate:treatment_group:biological_replicate) + (1 | Plate),
    data = df_HS6ST1
  )

# Plot

plot_HS6ST1 <- report_plot(fit_GPC4) +
  ggtitle("HS6ST1")


plot_HS6ST1_combined <- report_plot_combined(fit_GPC4) +
  ggtitle("HS6ST1")


## Plot grid

panel1 <- plot_grid(
  plot_AGG,
  NULL,
  plot_GPC1,
  plot_GPC4,
  plot_NDST1,
  plot_NDST2,
  plot_HS2ST1,
  plot_HS6ST1,
  ncol = 2,
  align = "h",
  axis = "t"
)

panel3 <- plot_grid(
  plot_AGG_combined,
  NULL,
  plot_GPC1_combined,
  plot_GPC4_combined,
  plot_NDST1_combined,
  plot_NDST2_combined,
  plot_HS2ST1_combined,
  plot_HS6ST1_combined,
  ncol = 2,
  align = "h",
  axis = "t"
)

# Panel 2: markers

## Plot 1: TUBB3

df_TUBB3 <- filter(df_big,
                 Detector == "TUBB3")

# Fit

fit_TUBB3 <-
  mixed(
    dCt_nexp_log2 ~ treatment_group * passage + cell_line + (1 |
                                                               Plate:treatment_group:biological_replicate) + (1 | Plate),
    data = df_TUBB3
  )

# Plot

plot_TUBB3 <- report_plot(fit_TUBB3) +
  ggtitle("TUBB3")


plot_TUBB3_combined <- report_plot_combined(fit_TUBB3) +
  ggtitle("TUBB3")


## Plot 2: NES

df_NES <- filter(df_big,
                   Detector == "NES")

# Fit

fit_NES <-
  mixed(
    dCt_nexp_log2 ~ treatment_group * passage + cell_line + (1 |
                                                               Plate:treatment_group:biological_replicate) + (1 | Plate),
    data = df_NES
  )

# Plot

plot_NES <- report_plot(fit_NES) +
  ggtitle("NES")


plot_NES_combined <- report_plot_combined(fit_NES) +
  ggtitle("NES")

## Plot 3: VIM

df_VIM <- filter(df_big,
                   Detector == "VIM")

# Fit

fit_VIM <-
  mixed(
    dCt_nexp_log2 ~ treatment_group * passage + cell_line + (1 |
                                                               Plate:treatment_group:biological_replicate) + (1 | Plate),
    data = df_VIM
  )

# Plot

plot_VIM <- report_plot(fit_VIM) +
  ggtitle("VIM")


plot_VIM_combined <- report_plot_combined(fit_VIM) +
  ggtitle("VIM")



## Plot 4: ADIPOQ

df_ADIPOQ <- filter(df_big,
                   Detector == "ADIPOQ")

# Fit

fit_ADIPOQ <-
  mixed(
    dCt_nexp_log2 ~ treatment_group * passage + cell_line + (1 |
                                                               Plate:treatment_group:biological_replicate) + (1 | Plate),
    data = df_ADIPOQ
  )

# Plot

plot_ADIPOQ <- report_plot(fit_ADIPOQ) +
  ggtitle("ADIPOQ")


plot_ADIPOQ_combined <- report_plot_combined(fit_ADIPOQ) +
  ggtitle("ADIPOQ")

## Plot grid

panel2 <- plot_grid(
  plot_TUBB3,
  plot_NES,
  plot_VIM,
  plot_ADIPOQ,
  ncol = 1,
  align = "h",
  axis = "t"
)

panel4 <- plot_grid(
  plot_TUBB3_combined,
  plot_NES_combined,
  plot_VIM_combined,
  plot_ADIPOQ_combined,
  ncol = 1,
  align = "h",
  axis = "t"
)

# Export plots

ggsave(filename = "./plots/MidYear_panel1.png",
       plot = panel1,
       scale = 1.1,
       width = 8,
       height = 10)

ggsave(filename = "./plots/MidYear_panel2.png",
       plot = panel2,
       scale = 1.1,
       width = 4,
       height = 10)

ggsave(filename = "./plots/MidYear_panel3.png",
       plot = panel3,
       scale = 1.1,
       width = 8,
       height = 10)

ggsave(filename = "./plots/MidYear_panel4.png",
       plot = panel4,
       scale = 1.1,
       width = 4,
       height = 10)


