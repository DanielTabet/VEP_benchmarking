library(data.table)
library(stringr)
library(tidyr)

# Load all phenotypes
measurements = fread('common/VEP_benchmarking_phenotypes.csv')

# Load config file
#source("lib/loadConfig.R")
#config = loadConfig(args[1])

# Set arguments
GENE = tolower(config$gene)
OUTPUT_PATH = = config$output_path

# Load participant eids
eids = readLines(sprintf(paste0(OUTPUT_PATH, "/%s/%s_unique_eids_filtered.txt"), toupper(GENE), GENE))

# Load phenotype selection
geneList = fread(config$gene_list)

bPhenotypes = geneList[gene_symbol == toupper(GENE),
                       .(field_code = str_split(phenotype_codes, "\n"),
                         field_description = str_split(phenotype_descriptions, "\n"))]

bPhenotypes = unnest(bPhenotypes, cols = c(field_code, field_description))
bPhenotypes = as.data.table(bPhenotypes)
bPhenotypes[, field_code := str_remove(field_code, "x")]
bPhenotypes[, field_code := str_remove(field_code, "\r")]
bPhenotypes[, field_id := as.numeric(str_split(field_code, "_", simplify = T)[, 1])]

# Attach phenotype category information
phenoCategory = fread(config$phenotype_list)
bPhenotypes = merge(bPhenotypes, phenoCategory, by = "field_id")

# Select relevant phenotypes and eids to match variants
colsToKeep <- c("eid", bPhenotypes$field_code)
measurements = measurements[, ..colsToKeep]
measurements = measurements[eid %in% eids]

# Save
outputFile = sprintf(paste0(OUTPUT_PATH, "/%s/%s_phenotypes_filtered.csv"), toupper(GENE), GENE)
fwrite(measurements, outputFile)
