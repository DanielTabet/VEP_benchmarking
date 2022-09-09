source("lib/packageCheck.R")

# Check argument
args = commandArgs(trailingOnly = T)
if (length(args) != 3) {
  cat("Usage:\n")
  cat("  Rscript worker.R <configuration-file> <gene-name> <log-file>")
  quit()
}

# Load configuration file
source("lib/loadConfig.R")
config = loadConfig(args[1])

# Add gene name and ensembl IDs to the config
config$gene = tolower(args[2])
config$ensembl_id = toupper(args[3])
config$canonical_transcript_id = toupper(args[4])

# Start logging
if (is.na(args[5])) stop("Log file path not specified")
library(logr)
logFile = log_open(args[5], logdir = F)

# Helper function: detach all packages
detachAllPackages = function() {
  basic.packages.blank = c(    
    "stats",    
    "graphics",    
    "grDevices",    
    "utils",   
    "datasets",  
    "methods",    
    "base"    
  )    
  basic.packages = paste("package:", basic.packages.blank, sep = "")   
  package.list = search()[ifelse(unlist(gregexpr("package:", search())) == 1, TRUE, FALSE)]   
  package.list = setdiff(package.list, basic.packages)   
  if (length(package.list) > 0) {   
    for (package in package.list) {   
      detach(package, character.only = TRUE)   
    }   
  }    
}

# Set config environment and run source
setUpEnvAndRun = function(prompt, file) {
  env = new.env()
  env$config = config
  
  # Run script
  logr::log_print(prompt)
  suppressPackageStartupMessages(sys.source(file, envir = env, keep.source = F))
  
  # Cleanup packages
  detachAllPackages()
}

# Create output folder
GENE = tolower(config$gene)
outputPath = paste("output", toupper(GENE), sep = "/")
if (!dir.exists(outputPath)) dir.create(outputPath)

# Parse variants
setUpEnvAndRun("Parsing variants", "lib/parseVariants.R")

# Prepare burden-test associated phenotypes
library(data.table)
library(stringr)
geneList = fread(config$gene_list)
bPhenotypes = geneList[gene_symbol == toupper(GENE),
                       .(field_code = str_split(phenotype_codes, "\n"),
                         field_description = str_split(phenotype_descriptions, "\n"))]
library(tidyr)
bPhenotypes = unnest(bPhenotypes, cols = c(field_code, field_description))
bPhenotypes = as.data.table(bPhenotypes)
bPhenotypes[, field_code := str_remove(field_code, "x")]
bPhenotypes[, field_id := as.numeric(str_split(field_code, "_", simplify = T)[, 1])]

# Attach phenotype category information
phenoCategory = fread(config$phenotype_list)
bPhenotypes = merge(bPhenotypes, phenoCategory, by = "field_id")

# Parse phenotypes
setUpEnvAndRun("Parsing phenotypes", "lib/parsePhenotypes.R")

# Merge phenotypes
setUpEnvAndRun("Merging phenotypes and variants", "lib/mergePhenotypesAndVariants.R")

# Benchmark variant effect predictors using quantitative traits
setUpEnvAndRun("Benchmarking variant predictors using quantitative traits", "lib/benchmarkWithQuantTraits.R")

# Benchmark variant effect predictors using categorical traits
setUpEnvAndRun("Benchmarking variant predictors using categorical traits", "lib/benchmarkWithCatTraits.R")

logr::log_close()