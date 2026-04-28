# Script to reinstall all required packages for genetic diversity analysis
# Run this script to install all necessary packages

# List of required packages
required_packages <- c(
  "here",
  "tidyverse",
  "ggpubr",
  "terra",
  "sf",
  "gridExtra",
  "rpaleoclim",
  "spdep",
  "spatialreg",
  "MuMIn",
  "cowplot",
  "performance",
  "gstat",
  "sp",
  "furrr",
  "interactions",
  "emmeans",
  "MetBrewer",
  "nlme",
  "broom.mixed",
  "emmeans",
  "jsonlite"  # Adding jsonlite since it was causing the error
)

# Function to install packages if not already installed
install_if_missing <- function(package_name) {
  if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
    cat("Installing package:", package_name, "\n")
    install.packages(package_name, dependencies = TRUE)
    library(package_name, character.only = TRUE)
  } else {
    cat("Package", package_name, "is already installed\n")
  }
}

# Install all packages
cat("Starting package installation...\n")
for (pkg in required_packages) {
  install_if_missing(pkg)
}

# Verify all packages are installed
cat("\nVerifying package installation...\n")
for (pkg in required_packages) {
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("✓", pkg, "is installed and loaded\n")
  } else {
    cat("✗", pkg, "failed to install\n")
  }
}

cat("\nPackage installation complete!\n")
