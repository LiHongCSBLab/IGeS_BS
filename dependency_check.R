# List of CRAN and Bioconductor packages
cran_packages <- c(
  "ggplot2", "patchwork", "dplyr", "readxl", "plot3D", "fmsb", "ggrepel",
  "fs", "plotrix", "plyr", "philentropy", "reshape2", "ggplotify", "grid",
  "cowplot", "parallel", "stringr", "data.table", "getopt"
)

bioc_packages <- c("waddR", "prada", "rhdf5", "limma")

# Install CRAN packages if missing
install_if_missing <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
  }
}

# Install Bioconductor packages if missing
install_bioc_if_missing <- function(pkgs) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg, ask = FALSE)
    }
  }
}

# Run installations
install_if_missing(cran_packages)
install_bioc_if_missing(bioc_packages)

