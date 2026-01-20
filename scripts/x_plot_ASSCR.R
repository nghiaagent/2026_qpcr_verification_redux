quant_limma <- readRDS(
    file = here(
        "output",
        "quant_limma.RDS"
    )
)

quant_big <- readRDS(
    file = here(
        "output",
        "quant_big.RDS"
    )
)

quant_small <- quant_big %>%
    filter(
        gene == "JUNB",
        treatment == "Control"
    )

ggplot(
    quant_small,
    aes(
        x = passage,
        y = dct_nexp_log2_offset,
        color = passage
    )
) +
    geom_boxplot() +
    geom_jitter(width = 0.1) +
    theme_classic() +
    labs(
        x = "Passage",
        y = "Relative log-expression"
    ) +
    scale_color_manual(values = c("#156082", "#196b24"))

ggsave(filename = "temp.png", width = 8, height = 6, scale = 0.6)
