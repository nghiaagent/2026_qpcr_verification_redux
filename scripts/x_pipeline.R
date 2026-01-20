# Run all scripts of pipeline

source(file.path(
    "scripts",
    "0_env_prep.R"
))

source(here(
    "scripts",
    "1_import.R"
))

source(here(
    "scripts",
    "2_clean.R"
))

source(here(
    "scripts",
    "3_load_limma.R"
))

source(here(
    "scripts",
    "4_analyse_limma.R"
))

source(here(
    "scripts",
    "5_plot_all.R"
))
