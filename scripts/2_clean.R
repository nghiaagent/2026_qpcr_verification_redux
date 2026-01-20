# Import packages
library("tidyverse")

# Load raw data
## Sort where appropriate
quant_qpcr <- readRDS(file = here("output", "quant_qpcr.RDS")) %>%
  arrange(
    sample_id,
    gene,
    biological_replicate,
    technical_replicate
  )

## Special filter: Remove GADD45B for Plate 5 Attempt 1
quant_qpcr <- quant_qpcr %>%
  filter(
    !(filename == "2024-11-26 MN TXome Verification Plate 5 Attempt 1" &
      gene == "GADD45B")
  )

# Set default values
## Threshold for qPCR filtering
threshold <- 2

## Plot specs
plot_specs <- list(
  scale = .9,
  width = 12,
  height = 8,
  units = "in",
  dpi = 300
)

# Clean 18S endogenous control
## Filter 18S values
quant_18S <- quant_qpcr %>%
  ### Select only 18S from dataset
  filter(gene == "18S") %>%
  ### Clean using custom logic to consider values more than 2 Ct out as outliers
  qPCR_clean_strict(
    Ct = ct,
    threshold = threshold,
    filename,
    sample_id,
    biological_replicate
  ) %>%
  ungroup() %>%
  ### Sort table for easy viewing
  arrange(
    filename,
    sample_id,
    passage,
    treatment,
    gene,
    technical_replicate
  )

## Average tech reps to get 18S value for each biol rep
quant_18S_summary <- quant_18S %>%
  group_by(name, filename) %>%
  summarise(ct_18S = mean(ct))

# Clean target genes
## Append 18S to target genes data
quant_qpcr_targets <- quant_qpcr %>%
  filter(!gene == "18S") %>%
  left_join(
    quant_18S_summary,
    by = join_by(
      name == name,
      filename == filename
    )
  )

## Clean data
## Same method as 18S
quant_qpcr_targets_clean <- quant_qpcr_targets %>%
  qPCR_clean_strict(
    Ct = ct,
    threshold = threshold,
    gene,
    name
  )

## View conditions with outliers present
quant_outliers <- anti_join(
  quant_qpcr_targets,
  quant_qpcr_targets_clean
)

# Find dCt
quant_big <- quant_qpcr_targets_clean %>%
  mutate(dct = ct - ct_18S) %>%
  mutate(dct_neg = max(dct) - dct) %>%
  mutate(dct_nexp = 2^-dct) %>%
  mutate(dct_nexp_log2 = log2(2^-dct)) %>%
  mutate(dct_nexp_log2_offset = log2(2^-dct) + 50) %>%
  ungroup()

# QC plots
## Ct histogram of all genes
quant_qpcr %>%
  mutate(
    ct_neg = (ct * -1) %>%
      replace_na(-60)
  ) %>%
  ggplot(
    aes(
      x = ct_neg,
      color = gene,
      fill = gene
    )
  ) +
  geom_density() +
  facet_grid(gene ~ cell_line)

ggsave(
  filename = here(
    "output",
    "plots_qc",
    "histogram_ct.png"
  ),
  scale = plot_specs$scale,
  width = plot_specs$width,
  height = plot_specs$height,
  units = plot_specs$units,
  dpi = plot_specs$dpi
)

## dCt histogram of all genes
quant_big %>%
  mutate(
    dct_neg = dct_neg %>%
      replace_na(-60)
  ) %>%
  ggplot(aes(
    x = dct_neg,
    color = gene,
    fill = gene
  )) +
  geom_density() +
  facet_grid(gene ~ cell_line)

ggsave(
  filename = here(
    "output",
    "plots_qc",
    "histogram_dct.png"
  ),
  scale = plot_specs$scale,
  width = plot_specs$width,
  height = plot_specs$height,
  units = plot_specs$units,
  dpi = plot_specs$dpi
)

# Export data
saveRDS(
  quant_big,
  file = here("output", "quant_big.RDS")
)
