################################################################################
############### Installation of required libriaries ############################
################################################################################
# Ref : Cornet et al. in prep
# Date: 2025-08-27
# Purpose:
#   Ensure that all required R packages are installed.
#   - Does NOT load the packages.
#   - Installs only those missing from the local environment.
#
################################################################################

# List of required packages
packages <- c(
  "bayesplot", "corrplot", "dendextend", "dplyr", "EnvStats",
  "ggnewscale","ggpattern", "ggplot2", "ggpubr", "ggrepel", "ggspatial", 
  "grid", "gridExtra", "locfit", "MASS", "patchwork", "png", "progress", 
  "purrr", "raster", "tidyr", "tidyverse"
)

# Identify packages that are not yet installed
to_install <- packages[!(packages %in% installed.packages()[,"Package"])]

# Install missing packages (if any)
if(length(to_install) > 0){
  install.packages(to_install)
}

if (!"cmdstanr" %in% installed.packages()[,"Package"]){
  install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
  cmdstanr::install_cmdstan()
}
