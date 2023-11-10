CRAN_LIBRARIES = c("yaml", "data.table", "dplyr", "tidyr", "stringr",
                   "ggtext", "glue", "forcats", "gridExtra", "fs", "Rmisc",
                   "reshape2", "stringr", "logr", "BiocManager", "qvalue", "remotes", "retry")
#BIOCONDUCTOR_LIBRARIES = c("EnsDb.Hsapiens.v86")
GITHUB_LIBRARIES = c("jweile/yogiroc")

# Install needed libraries
cat("Checking missing libraries\n")
missingLibs = CRAN_LIBRARIES[!CRAN_LIBRARIES %in% installed.packages()[, "Package"]]
if (length(missingLibs)) {
  cat("Installing missing CRAN libraries:", paste0(missingLibs, collapse = ", "), "\n")
  install.packages(missingLibs, lib = Sys.getenv("R_LIBS_USER"))
}

missingLibs = BIOCONDUCTOR_LIBRARIES[!BIOCONDUCTOR_LIBRARIES %in% installed.packages()[, "Package"]]
if (length(missingLibs)) {
  cat("Installing missing Bioconductor libraries:", paste0(missingLibs, collapse = ", "), "\n")
  BiocManager::install(missingLibs, lib = Sys.getenv("R_LIBS_USER"))
}

library(stringr)
missingLibs = str_split(GITHUB_LIBRARIES, fixed("/"), simplify = T)[,2]
missingLibs = GITHUB_LIBRARIES[which(!missingLibs %in% installed.packages()[, "Package"])]
if (length(missingLibs)) {
  cat("Installing missing Github libraries:", paste0(missingLibs, collapse = ", "), "\n")
  remotes::install_github(missingLibs, lib = Sys.getenv("R_LIBS_USER"))
}