# Import packages
library("limma")
library("tidyverse")

# Load data
quant_limma <- readRDS(
  file = here(
    "output",
    "quant_limma.RDS"
  )
)

# Create design matrix that examines effect of condition
## Create design matrix / define model
design <- model.matrix(~sample_id, quant_limma$targets)
colnames(design) <- colnames(design) %>%
  make.names()

## Define contrasts
matrix_contrasts <- makeContrasts(
  # Coefs 1 - 5: hMSC coefficients
  ## Effect of passage
  P13vsP5_Ctr = sample_idhMSC.21558.P13.D5.Control - 0,
  P13vsP5_Hep = sample_idhMSC.21558.P13.D5.Heparin -
    sample_idhMSC.21558.P5.D5.Heparin,

  ## Effect of heparin
  Hep_hMSC_P5 = sample_idhMSC.21558.P5.D5.Heparin - 0,
  Hep_hMSC_P13 = sample_idhMSC.21558.P13.D5.Heparin -
    sample_idhMSC.21558.P13.D5.Control,

  ## Effect of passage x heparin
  Hep_x_Passage_hMSC = (sample_idhMSC.21558.P13.D5.Heparin -
    sample_idhMSC.21558.P13.D5.Control) -
    (sample_idhMSC.21558.P5.D5.Heparin - 0),

  # Coefs 6 - 9: Effect of Hep in other cell pops
  ## Coef 6: RCX
  Hep_RCX = sample_idReNcell.CX.D5.Heparin - sample_idReNcell.CX.D5.Control,

  ## Coef 7: RVM
  Hep_RVM = sample_idReNcell.VM.D5.Heparin - sample_idReNcell.VM.D5.Control,

  ## Coef 8: 1321N1
  Hep_1321N1 = sample_id1321N1.D5.Heparin - sample_id1321N1.D5.Control,

  ## Coef 9: SH
  Hep_SH = sample_idSH.SY5Y.D3.Heparin - sample_idSH.SY5Y.D3.Control,

  # Coef 10 - 19: Compare between cell pops (baseline)
  ## Other vs. hMSC
  ### RCX vs. hMSC
  RCX_vs_hMSC = sample_idReNcell.CX.D5.Control - 0,
  ### RVM vs. hMSC
  RVM_vs_hMSC = sample_idReNcell.VM.D5.Control - 0,
  ### 1321 vs. hMSC
  X1321_vs_hMSC = sample_id1321N1.D5.Control - 0,
  ### SH-SY5Y vs. hMSC
  SH_vs_hMSC = sample_idSH.SY5Y.D3.Control - 0,

  ## Other vs. RCX
  ### RVM vs. RCX
  RVM_vs_RCX = sample_idReNcell.VM.D5.Control - sample_idReNcell.CX.D5.Control,
  ### 1321N1 vs. RCX
  X1321_vs_RCX = sample_id1321N1.D5.Control - sample_idReNcell.CX.D5.Control,
  ### SH-SY5Y vs. RCX
  SH_vs_RCX = sample_idSH.SY5Y.D3.Control - sample_idReNcell.CX.D5.Control,

  ## Other vs. RVM
  ### 1321 vs. RVM
  X1321_vs_RVM = sample_id1321N1.D5.Control - sample_idReNcell.VM.D5.Control,
  ### SH-SY5Y vs. RVM
  SH_vs_RVM = sample_idSH.SY5Y.D3.Control - sample_idReNcell.VM.D5.Control,

  ## Other vs. 1321N1
  ### SH-SY5Y vs. 1321N1
  SH_vs_X1321 = sample_idSH.SY5Y.D3.Control - sample_id1321N1.D5.Control,

  # Coef 20 - 29: Compare between cell pops (response to heparin)
  ## Other vs. hMSC
  ### RCX vs. hMSC
  RCX_vs_hMSC_treat = (sample_idReNcell.CX.D5.Heparin -
    sample_idReNcell.CX.D5.Control) -
    (sample_idhMSC.21558.P5.D5.Heparin - 0),
  ### RVM vs. hMSC
  RVM_vs_hMSC_treat = (sample_idReNcell.VM.D5.Heparin -
    sample_idReNcell.VM.D5.Control) -
    (sample_idhMSC.21558.P5.D5.Heparin - 0),
  ### 1321 vs. hMSC
  X1321_vs_hMSC_treat = (sample_id1321N1.D5.Heparin -
    sample_id1321N1.D5.Control) -
    (sample_idhMSC.21558.P5.D5.Heparin - 0),
  ### SH-SY5Y vs. hMSC
  SH_vs_hMSC_treat = (sample_idSH.SY5Y.D3.Heparin -
    sample_idSH.SY5Y.D3.Control) -
    (sample_idhMSC.21558.P5.D5.Heparin - 0),

  ## Other vs. RCX
  ### RVM vs. RCX
  RVM_vs_RCX_treat = (sample_idReNcell.VM.D5.Heparin -
    sample_idReNcell.VM.D5.Control) -
    (sample_idReNcell.CX.D5.Heparin - sample_idReNcell.CX.D5.Control),
  ### 1321N1 vs. RCX
  X1321_vs_RCX_treat = (sample_id1321N1.D5.Heparin -
    sample_id1321N1.D5.Control) -
    (sample_idReNcell.CX.D5.Heparin - sample_idReNcell.CX.D5.Control),
  ### SH-SY5Y vs. RCX
  SH_vs_RCX_treat = (sample_idSH.SY5Y.D3.Heparin -
    sample_idSH.SY5Y.D3.Control) -
    (sample_idReNcell.CX.D5.Heparin - sample_idReNcell.CX.D5.Control),

  ## Other vs. RVM
  ### 1321 vs. RVM
  X1321_vs_RVM_treat = (sample_id1321N1.D5.Heparin -
    sample_id1321N1.D5.Control) -
    (sample_idReNcell.VM.D5.Heparin - sample_idReNcell.VM.D5.Control),
  ### SH-SY5Y vs. RVM
  SH_vs_RVM_treat = (sample_idSH.SY5Y.D3.Heparin -
    sample_idSH.SY5Y.D3.Control) -
    (sample_idReNcell.VM.D5.Heparin - sample_idReNcell.VM.D5.Control),

  ## Other vs. 1321
  ### SH vs. 1321
  SH_vs_X1321_treat = (sample_idSH.SY5Y.D3.Heparin -
    sample_idSH.SY5Y.D3.Control) -
    (sample_id1321N1.D5.Heparin - sample_id1321N1.D5.Control),
  levels = design
)

rownames(matrix_contrasts) <- colnames(design)

# Fit with limma
## Fit model
fit <- quant_limma %>%
  lmFit(design) %>%
  eBayes()

fit_contrasts <- fit %>%
  contrasts.fit(matrix_contrasts) %>%
  eBayes()

# Show results
## Test summary
summary(decideTests(fit_contrasts))

# Save data
saveRDS(
  fit_contrasts,
  file = here(
    "output",
    "fit_contrasts.RDS"
  )
)
