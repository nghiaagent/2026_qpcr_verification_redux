library("colorspace")

# Define colour palettes

# Palette 1: Set of colours for cell lines and timepoints
# hMSC P5, hMSC P7, hMSC P13, RCX, RVM, 1321, SH (7 colours) + treatment variant
# For Figs 1C, 2D, 3G, 5
## Get colours for palettes
## Extract colours from R Okabe-Ito palette in desired order)
## Consistent with Volcano3D plot
### P5: Blue
### P7: Red (Vermillion)
### P13: Green
### RCX: Yellow
### RVM: Orange
### 1321: Purple
### SH: Grey
palette <- c(
    palette.colors()[[6]],
    palette.colors()[[7]],
    palette.colors()[[4]],
    palette.colors()[[5]],
    palette.colors()[[2]],
    palette.colors()[[8]],
    palette.colors()[[9]]
)

palette_merge <- vec_interleave(
    lighten(palette, 0.2),
    darken(palette, 0.2),
)
