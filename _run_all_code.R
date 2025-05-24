rm(list=ls()) # Clear work space
system("rm _Outputs/*")
system("rm _R_Data_Files/*")

# List of required packages
packages <- c(
  "betareg", "binom", "biostat3", "broom", "dplyr", "distributionsrd", "ggplot2",
  "ggthemes", "glmmTMB", "kableExtra", "knitr", "lmtest", "lubridate", "magrittr",
  "MASS", "patchwork", "pracma", "stringr", "tidyverse", "tidyr"
)

# Install any packages that are missing
installed <- packages %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(packages[!installed])
}
# Load all packages
invisible(lapply(packages, library, character.only = TRUE))

# Run all scripts
source("01_calc_pvalues.R")
source("02_do_fits.R")
source("03_make_fits_table.R")
source("04_prop_comparisons.R")
source("05_stat_analysis.R")
source("06_simulations.R")
source("07_plot_simulations.R")
source("08_make_fig2.R")
source("09_phackfit.R")
source("10_make_figure3.R")
source("11_revisions_nonlinear.R")
source("12_revisions_quantile_regression.R")
source("13_revisions_changing_power&phacking.R")

