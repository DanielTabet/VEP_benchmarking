library(data.table)
library(stringr)

setwd("../")

###
# Count genes and unique variants
###

# Load genes
geneInfo = fread("common/genes.csv")
genes = unique(geneInfo$gene_symbol)

# For each gene, count unique variants
varCount = sapply(genes, function(gene) {
  variantPath = sprintf("output/%s/%s_variants_filtered.csv", toupper(gene), gene)
  if (!file.exists(variantPath)) return(NA)
  
  # Load variants
  variants = fread(variantPath)
  return(length(unique(variants$variant_name)))
})
cat("Total genes:", length(na.omit(varCount)), "\n")
cat("Total unique variants:", sum(varCount, na.rm = T), "\n")

###
# Count unique traits
###

# Load traits
traits = str_split(geneInfo$phenotype_codes, "\n")
traits = unique(unlist(traits))
cat("Total unique traits:", length(traits), "\n")

###
# Count unique UKB participants
###

# For each gene, count UKB participants
partCount = sapply(genes, function(gene) {
  variantPath = sprintf("output/%s/%s_variants_filtered.csv", toupper(gene), gene)
  if (!file.exists(variantPath)) return(NA)
  
  # Load variants
  variants = fread(variantPath)
  return(unique(variants$eid))
})
partCount = unique(unlist(partCount))
cat("Total unique UKB participants:", length(partCount), "\n")

###
# Count computational predictor coverage
###
# For each gene, count the number of scores predicted
scoreCount = sapply(genes, function(gene, predictor) {
  variantPath = sprintf("output/%s/%s_variants_filtered.csv", toupper(gene), gene)
  if (!file.exists(variantPath)) return(NA)
  
  # Load variants
  variants = fread(variantPath)
  
  predictorPath = sprintf("input/%s/allv_withweights.csv", toupper(gene))
  if (!file.exists(predictorPath)) return(NA)
  
  # Load predictor scores
  predictors = fread(predictorPath)
  
  # Select filtered variants and predictor
  predictors = predictors[variant_name %in% variants$variant_name]
  predictors = predictors[, predictor, with = F]
  
  return(length(na.omit(unlist(predictors, use.names = F))))
}, predictor = "evm_score_flipped")
cat("The number of variants for which EVmutation made predictions:", sum(scoreCount, na.rm = T), "\n")
