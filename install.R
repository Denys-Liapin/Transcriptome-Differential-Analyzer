message("Starting package check and installation...")

# Install BiocManager if not present
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# List of CRAN packages
cran_packages <- c(
  "readxl",
  "ggplot2",
  "ggvenn",
  "ggpubr",
  "pheatmap",
  "RColorBrewer",
  "stringr",
  "ggrepel"
)

# Check and install CRAN packages
new_cran_packages <- cran_packages[!(cran_packages %in% installed.packages()[,"Package"])]
if (length(new_cran_packages) > 0) {
  message("Installing CRAN packages: ", paste(new_cran_packages, collapse = ", "))
  install.packages(new_cran_packages)
} else {
  message("All CRAN packages are already installed.")
}

# List of Bioconductor packages
bioc_packages <- c(
  "DESeq2",
  "EnhancedVolcano",
  "clusterProfiler",
  "org.Dm.eg.db",
  "KEGGREST"
)

# Check and install Bioconductor packages
new_bioc_packages <- bioc_packages[!(bioc_packages %in% installed.packages()[,"Package"])]
if (length(new_bioc_packages) > 0) {
  message("Installing Bioconductor packages: ", paste(new_bioc_packages, collapse = ", "))
  BiocManager::install(new_bioc_packages, update = FALSE)
} else {
  message("All Bioconductor packages are already installed.")
}

message("Done! All required packages have been successfully installed.")