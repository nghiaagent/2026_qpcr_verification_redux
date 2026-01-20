# Import packages
library("tidyverse")

# Import all qPCR data as tsv
## 1. Get sample list from sample table
## 2. Import files from folder, _MANUALLY_ remove header and footer from generated csv file.
## 3. Filter for samples with included = TRUE
## 4. Clean up dataset
## 5. Export

## NOTE: When exporting from QuantStudio,
## export only results tab,
## file format: .txt

# Set sample table path
path_sample <- here::here(
  "input",
  "annotation",
  "2024_hmsc_qpcr_verification_sample_table.csv"
)

# Construct sample table
table_samples <- read_csv(path_sample) %>%
  ## Recode biological variables as factors
  ## To allow ordering of levels of variables.
  mutate(
    cell_line = cell_line %>%
      factor(
        levels = c(
          "hMSC-21558",
          "ReNcellCX",
          "ReNcellVM",
          "1321N1",
          "SH-SY5Y"
        )
      ),
    passage = passage %>%
      factor(levels = c("P5", "P13")),
    day = day %>%
      factor(levels = c("D3", "D5")),
    treatment = treatment %>%
      factor(levels = c("Control", "Heparin"))
  ) %>%
  ## Sort table for easy viewing and recoding of more variables as factors
  arrange(
    cell_line,
    filename,
    passage,
    day,
    treatment
  ) %>%
  ## Recode condition and sample name as factors
  mutate(
    sample_id = sample_id %>%
      factor() %>%
      fct_inorder(),
    name = name %>%
      factor() %>%
      fct_inorder()
  ) %>%
  ## Group by plate (filename) for filtering of biological replicates
  arrange(filename) %>%
  group_by(filename)

# Define list of file  paths, biological replicates to include
## Establish list of biological replicates to be included from each file
samples_include <- table_samples %>%
  group_split() %>%
  setNames(group_keys(table_samples)[[1]]) %>%
  map(\(x) filter(x, included_in_dataset == TRUE)$name)

## List of file paths
list_files <- here::here(
  "input",
  "txt",
  str_c(table_samples$filename, ".txt")
) %>%
  unique() %>%
  setNames(unique(table_samples$filename))

## Check all files exist
map(
  list_files,
  \(file) {
    if (file.exists(file) == FALSE) {
      file.exists(file)
      stop(str_wrap(
        "There is a filename mismatch between
                sample sheet and data"
      ))
    } else {
      file.exists(file)
    }
  }
)

# Import files. Skip csv headers manually
quant_qpcr_raw <- list_files %>%
  imap(\(filepath, filename) {
    read_tsv(
      filepath,
      skip = 53
    ) %>%
      mutate(filename = filename)
  })

# Clean raw data
quant_qpcr <- quant_qpcr_raw %>%
  ## Filter for samples marked as TRUE on sample table in each file
  map2(
    .,
    samples_include,
    \(data, include) {
      filter(data, `Sample Name` %in% include)
    }
  ) %>%
  ## Merge list of dataframes into one
  bind_rows() %>%
  ## Select only relevant columns
  dplyr::select(
    `Sample Name`,
    `Target Name`,
    `CT`,
    `filename`
  ) %>%
  ## Append sample data
  right_join(
    table_samples,
    .,
    by = join_by(
      `name` == `Sample Name`,
      `filename` == `filename`
    )
  ) %>%
  ## Remove grouping based on filename applied earlier
  ungroup() %>%
  ## Convert CT to numeric
  mutate(
    CT = CT %>%
      as.numeric()
  ) %>%
  ## Add biological and tech rep numbers
  mutate(
    biological_replicate = biological_replicate %>%
      as.character(),
    technical_replicate = rep(1:4, nrow(.) / 4) %>%
      as.character() %>%
      str_c(
        biological_replicate,
        .,
        sep = "_"
      )
  ) %>%
  relocate(
    technical_replicate,
    .after = biological_replicate
  ) %>%
  ## Rename variables to fit style guide
  dplyr::rename(
    "gene" = "Target Name",
    "ct" = "CT"
  )


## Save data
saveRDS(
  quant_qpcr,
  file = here(
    "output",
    "quant_qpcr.RDS"
  )
)

saveRDS(
  table_samples,
  file = here(
    "output",
    "table_samples.RDS"
  )
)
