library(data.table)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggtext)

# Check argument
args = commandArgs(trailingOnly = T)
if (length(args) != 2) {
  cat("Usage:\n")
  cat("  Rscript tools/plotPRScoreDistri.R <configuration-file> <gene-name> <predictor>")
  quit()
}

# Load configuration file
source("lib/loadConfig.R")
config = loadConfig(args[1])

# Get gene name and field id
GENE = tolower(args[2])
PREDICTOR = args[3]

# Load phenotype descriptions
phenoDescription = fread(config$phenotype_list)
geneList = fread(config$gene_list)
bPhenotypes = geneList[gene_symbol == toupper(GENE),
                       .(field_code = str_split(phenotype_codes, "\n"),
                         field_description = str_split(phenotype_descriptions, "\n"))]
bPhenotypes = unnest(bPhenotypes, cols = c(field_code, field_description))
bPhenotypes = as.data.table(bPhenotypes)
bPhenotypes[, field_code := str_remove(field_code, "x")]
bPhenotypes[, field_id := as.numeric(str_split(field_code, "_", simplify = T)[, 1])]
bPhenotypes = merge(bPhenotypes, phenoDescription, by = "field_id")

# Plot variant score distributions
phenotypes = bPhenotypes[field_type == "categorical", as.character(field_id)]
plots = lapply(1:length(phenotypes), function(index, predictor) {
  # Get phenotype id
  phenotype = phenotypes[[index]]
  
  # Load precision recall analysis output
  filePath = sprintf("output/%s/%s_%s_auprc.rds", toupper(GENE), GENE, phenotype)
  if (!file.exists(filePath))
    stop(sprintf("Missing Precision Recall analysis output: %s", filePath))
  results = readRDS(filePath)
  
  # Load variants
  vTable = results$variants
  
  # Plot distribution
  plotTable = vTable[, .(p = !(is.na(p) | p == ""), score = get(predictor))]
  plotTable[score > 1, score := 1]
  cuts = cut_interval(plotTable$score, 11)
  breaks = seq(round(min(plotTable$score), 0), round(max(plotTable$score), 0), length.out = 11)
  breaks = as.character(round(as.numeric(breaks), 1))
  names(breaks) = levels(cuts)
  breaks[11] = ">=1"
  plotTable$bin = as.numeric(cuts)
  plotTable = plotTable[, .(perc = as.double(.N)), by = c("p", "bin")]
  plotTable[p == TRUE, c("n", "perc") := list(sum(perc), perc / sum(perc))]
  plotTable[p == FALSE, c("n", "perc") := list(sum(perc), perc / sum(perc))]
  plotTable[, se := sqrt(perc * (1 - perc)/n)]
  
  plot = ggplot(plotTable, aes(x = bin, y = perc, color = p)) +
    geom_point(aes(shape = p), size = 2) + geom_line(size = 1) +
    geom_errorbar(aes(ymin = perc - se, ymax = perc + se), width = 0.2) +
    scale_x_discrete(breaks = names(breaks), limits = names(breaks), labels = breaks) +
    scale_y_continuous(labels = scales::percent) +
    scale_color_manual(values = c("#457b9d", "#e76f51")) +
    labs(x = "Participant's sum of commonly-scaled variant scores", y = "% Participants",
         color = bPhenotypes[field_id == phenotype, field_description],
         shape = bPhenotypes[field_id == phenotype, field_description]) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom", legend.title = element_markdown(),
          axis.text = element_text(size = 12))
  
  return(plot)
}, predictor = PREDICTOR)
plots = arrangeGrob(grobs = plots,
                    top = grid::textGrob(toupper(GENE),
                                         gp = grid::gpar(fontsize = 20)))
filePath = sprintf("output/%s/%s_%s_variant_scores_dist.png", toupper(GENE), GENE, PREDICTOR)
ggsave(filePath, plots, width = 18, height = 18, units = "in")
