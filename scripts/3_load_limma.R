# Import packages
library("tidyverse")

# Load data
table_samples <- readRDS(
  file = here(
    "output",
    "table_samples.RDS"
  )
)

quant_big <- readRDS(
  file = here(
    "output",
    "quant_big.RDS"
  )
)

# Load into limma
## Convert df into matrix
df_matrix <- quant_big %>%
  mutate(
    col_name = str_c(
      sample_id,
      technical_replicate,
      sep = "_"
    )
  ) %>%
  dplyr::select(gene, dct_nexp_log2, col_name) %>%
  pivot_wider(
    names_from = col_name,
    values_from = dct_nexp_log2
  )

quant_limma <- as.matrix(df_matrix[-1])
row.names(quant_limma) <- as_vector(df_matrix[1])

## Generate design matrix for limma
### Create targets table
targets <- quant_limma %>%
  colnames() %>%
  as_tibble() %>%
  mutate(sample = value) %>%
  separate(
    value,
    into = c(
      "sample_id",
      "biological_replicate",
      "technical_replicate"
    ),
    sep = "_"
  ) %>%
  left_join(
    table_samples %>%
      ungroup() %>%
      select(
        sample_id,
        cell_line,
        passage,
        day,
        treatment
      ) %>%
      unique(),
    by = join_by(
      "sample_id" == "sample_id"
    ),
    relationship = "many-to-one"
  ) %>%
  mutate(
    sample_id = sample_id %>%
      factor() %>%
      fct_inorder(),
    batch = cell_line %>%
      case_match(
        "hMSC-21558" ~ "4",
        "ReNcellVM" ~ "5",
        "ReNcellCX" ~ "6",
        "SH-SY5Y" ~ "5",
        "1321N1" ~ "5"
      ) %>%
      factor() %>%
      fct_inorder()
  )

## Create EList object
dge <- new("EList")
dge$E <- quant_limma
dge$targets <- targets

### Add annotation data
dge$genes <- rownames(dge) %>%
  as_tibble() %>%
  dplyr::rename(GENENAME = value)

# Save data
saveRDS(
  dge,
  file = here(
    "output",
    "quant_limma.RDS"
  )
)
