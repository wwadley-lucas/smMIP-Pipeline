# CHIP Analysis Dependencies Installer
# Run this script to install all required packages

cat("Installing CHIP Analysis dependencies...\n")

# CRAN packages
cran_packages <- c(
  "optparse",       # Command line argument parsing
  "tidyverse",      # Data manipulation
  "rmarkdown",      # Report generation
  "knitr",          # R Markdown support
  "kableExtra",     # Table formatting
  "rentrez",        # NCBI/ClinVar API
  "httr",           # HTTP requests
  "jsonlite",       # JSON parsing
  "circlize",       # Circular visualizations
  "RColorBrewer",   # Color palettes
  "scales"          # Scale functions for ggplot2
)

for (pkg in cran_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, repos = "https://cloud.r-project.org")
  } else {
    cat(pkg, "already installed.\n")
  }
}

# BiocManager for Bioconductor packages
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# Bioconductor packages
bioc_packages <- c("ComplexHeatmap")

for (pkg in bioc_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Installing", pkg, "from Bioconductor...\n")
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  } else {
    cat(pkg, "already installed.\n")
  }
}

# TinyTeX for PDF generation
if (!require("tinytex", quietly = TRUE)) {
  install.packages("tinytex", repos = "https://cloud.r-project.org")
}

# Check if TinyTeX is installed
if (!tinytex::is_tinytex()) {
  cat("\nTinyTeX not found. Installing TinyTeX for PDF generation...\n")
  cat("This may take several minutes.\n")
  tinytex::install_tinytex()
} else {
  cat("TinyTeX already installed.\n")
}

cat("\nAll dependencies installed successfully!\n")
cat("You can now run: source('generate_patient_reports.R')\n")
