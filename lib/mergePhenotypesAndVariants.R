library(data.table)

# Set arguments if not passed in already as config
GENE = tolower(config$gene)
OUTPUT_PATH = = config$output_path

# Load filtered phenotypes and variants
phenotypes = fread(sprintf(paste0(OUTPUT_PATH, "/%s/%s_phenotypes_filtered.csv"), toupper(GENE), GENE))
variants = fread(sprintf(paste0(OUTPUT_PATH, "/%s/%s_variants_filtered.csv"), toupper(GENE), GENE))

# Merge variants and phenotypes
variants$eid = as.numeric(variants$eid)
variants = variants[!is.na(eid)]
variants = merge(variants, phenotypes, by = "eid")

# Select quantitative phenotypes
columns = colnames(variants)
columns = columns[!columns %in% bPhenotypes[field_type == "categorical", field_id]]
quantPhenotype = variants[, ..columns]
quantPhenotype[quantPhenotype == ""] = NA

# Save quantitative and all phenotypes separately
fwrite(quantPhenotype, sprintf(paste0(OUTPUT_PATH, "/%s/%s_quant-variants_measures.csv"), toupper(GENE), GENE))
fwrite(variants, sprintf(paste0(OUTPUT_PATH, "/%s/%s_variants_phenotypes.csv"), toupper(GENE), GENE))