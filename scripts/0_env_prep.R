# Load Windows fonts for use in project
library("extrafont")
loadfonts(device = "win")

# Define functions
library("here")

source(here::here(
    "scripts",
    "0_define_qPCR_outlier_detect.R"
))

source(here::here(
    "scripts",
    "0_define_results_clip.R"
))

source(here::here(
    "scripts",
    "0_define_colours.R"
))

# Define genes of interest
mndata_genes <- c(
    "GADD45B",
    "JUNB",
    "MAP3K7",
    "NFKB1"
)
