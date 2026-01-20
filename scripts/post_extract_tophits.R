# Load data

source("./scripts/analyse_limma.R")

# Extract logFC - treatment

fit_small <-
  fit_contrasts[fit_contrasts$genes$GENENAME %in% genes_GOI,]

topTable <- list(
  "P5_D3" = topTable(
    fit_small,
    coef = 1,
    number = Inf,
    sort.by = "none"
  ),
  "P5_D5" = topTable(
    fit_small,
    coef = 4,
    number = Inf,
    sort.by = "none"
  ),
  "P7_D3" = topTable(
    fit_small,
    coef = 2,
    number = Inf,
    sort.by = "none"
  ),
  "P7_D5" = topTable(
    fit_small,
    coef = 5,
    number = Inf,
    sort.by = "none"
  ),
  "P13_D3" = topTable(
    fit_small,
    coef = 3,
    number = Inf,
    sort.by = "none"
  ),
  "P13_D5" = topTable(
    fit_small,
    coef = 6,
    number = Inf,
    sort.by = "none"
  )
)

write.xlsx(topTable,
           file = "topTable.xlsx",
           asTable = TRUE)
