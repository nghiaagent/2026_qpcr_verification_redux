## Define function to clip logFC and padj in limma topTable to desired point
## For plotting with EnhancedVolcano
## Add relevant metadata
clip_results <- function(
    toptable,
    cutoff_logfc = 2.5,
    cutoff_padj = 1e-20,
    alpha = 0.05) {
    cutoff_logfc_neg <- cutoff_logfc * -1

    toptable <- tibble(toptable)

    # Clip logFC to the threshold
    toptable$logFC <- toptable$logFC %>%
        case_when(
            . >= cutoff_logfc ~ cutoff_logfc,
            . <= cutoff_logfc_neg ~ cutoff_logfc_neg,
            .default = .
        )

    # Clip padj to the threshold
    toptable$adj.P.Val <- toptable$adj.P.Val %>%
        case_when(
            . <= cutoff_padj ~ cutoff_padj,
            .default = .
        )

    # Add custom shapes to dots to identify clipped genes
    ## Normal dots: Shape 19
    ## Clipped (positive): Shape -9658
    ## Clipped (negative): Shape -9668
    ## Add names for legend
    toptable$volcano_shape <- case_when(
        toptable$logFC >= cutoff_logfc ~ -9658,
        toptable$logFC <= cutoff_logfc_neg ~ -9668,
        toptable$adj.P.Val <= cutoff_padj ~ 17,
        .default = 19
    )

    names(toptable$volcano_shape) <- case_when(
        toptable$volcano_shape == -9658 ~ str_c("logFC >", cutoff_logfc),
        toptable$volcano_shape == -9668 ~ str_c("logFC < -", cutoff_logfc),
        toptable$volcano_shape == 17 ~ str_c("padj <", cutoff_padj),
        toptable$volcano_shape == 19 ~ "Unclipped"
    )

    # Make clipped genes larger
    toptable$volcano_size <- case_when(
        toptable$logFC >= cutoff_logfc ~ 3,
        toptable$logFC <= cutoff_logfc_neg ~ 3,
        toptable$adj.P.Val <= cutoff_padj ~ 3,
        .default = 1
    )

    # Create new AveExpr column for MA plot
    toptable$AveExpr_new <- 1 / (toptable$AveExpr + 1)

    # Create new significance status level column
    toptable$colour <- case_when(
        toptable$adj.P.Val <= alpha ~ "red2",
        .default = "grey30"
    )

    names(toptable$colour) <- case_when(
        toptable$colour == "red2" ~ str_c("padj <", alpha),
        toptable$colour == "grey30" ~ "ns"
    )

    # Return modified data
    return(toptable)
}
