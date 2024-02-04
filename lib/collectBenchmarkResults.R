library(data.table)
library(stringr)
library(ggplot2)
library(forcats)
library(tidyr)
library(ggtext)
library(gridExtra)
library(qvalue)

# Load configuration file
source("lib/loadConfig.R")
config = loadConfig("config.yaml")

# Set arguments if not passed in already as config
VARIANT_PREDICTORS = config$variant_predictors

# Load genes
geneList = fread(config$gene_list)
genes = geneList$gene_symbol

# Helper function: extract AUPRC to predictor label
extractAUPRC = function(labels) {
  return(as.numeric(str_extract_all(labels, "\\d+\\.\\d+")))
}

# Load phenotype categories
phenoCategory = fread(config$phenotype_list)

# Calculate P values
pVals = lapply(1:length(genes), function(index) {
  gene = tolower(genes[index])
  cat(sprintf("\rProcess gene %s (%d/%d)...", gene, index, length(genes)))
  
  # Load phenotypes
  bPhenotypes = geneList[gene_symbol == toupper(gene),
                         .(field_code = str_split(phenotype_codes, "\n"),
                           field_description = str_split(phenotype_descriptions, "\n"))]
  bPhenotypes = unnest(bPhenotypes, cols = c(field_code, field_description))
  bPhenotypes = as.data.table(bPhenotypes)
  bPhenotypes[, field_code := str_remove(field_code, "x")]
  bPhenotypes[, field_code := str_remove(field_code, "\r")]
  bPhenotypes[, field_id := as.numeric(str_split(field_code, "_", simplify = T)[, 1])]
  
  # Attach phenotype category information
  bPhenotypes = merge(bPhenotypes, phenoCategory, by = "field_id")
  
  # Init p value variables
  pccPvals = NULL
  prcPvals = NULL
  
  # Check if correlations exist
  filePath = sprintf(paste0(OUTPUT_PATH, "/%s/%s_correlations.csv"), toupper(gene), gene)
  
  if (file.exists(filePath)) {
    # Load correlations
    corrs = fread(filePath)
    
    # Check if pccs exist
    if (!("pcc" %in% colnames(corrs))) return(NULL)
    
    # Remove NA columns
    corrs = corrs[complete.cases(corrs)]
    if (nrow(corrs) < 1) return(NULL)
    
    # Compute P values for each phenotype
    phenotypes = unique(corrs$field_id)
    pccPvals = lapply(phenotypes, function(phenotype) {
      # Compute P values for each predictor pair
      pccPvals = expand.grid(pred1 = corrs$name, pred2 = corrs$name)
      pccPvals = as.data.table(pccPvals)
      pccPvals$p_val = apply(pccPvals, 1, function(row) {
        # Get predictors
        pred1 = row[["pred1"]]
        pred2 = row[["pred2"]]
        
        # Get PCC distributions
        pred1Dist = as.numeric(str_split(corrs[name == pred1 & field_id == phenotype, pcc], fixed("|"), simplify = T)[1,])
        pred2Dist = as.numeric(str_split(corrs[name == pred2 & field_id == phenotype, pcc], fixed("|"), simplify = T)[1,])
        
        # Compute difference distribution
        diffDist = pred1Dist^2 - pred2Dist^2
        
        # Calculate the probability of P1 is not bigger than P2
        # (i.e. the difference is smaller than or equal to 0)
        pVal = sum(diffDist <= 0) / length(diffDist)
        
        return(pVal)
      })
      
      # Helper function: attach PCC to predictor name
      attachPCC = function(formName, fieldId) {
        val = corrs[name == formName & field_id == fieldId, avg_pcc]
        return(sprintf("%s (%.3f)", formName, val))
      }
      
      # Format labels
      pccPvals$pred1 = sapply(pccPvals$pred1, attachPCC, phenotype)
      pccPvals$pred2 = sapply(pccPvals$pred2, attachPCC, phenotype)
      
      # Attach information and return
      pccPvals$code = phenotype
      pccPvals$type = "PCC"
      return(pccPvals)
    })
    pccPvals = rbindlist(pccPvals)
  }
  
  # Check if AUPRC scores exist
  catePhenotypes = bPhenotypes[field_type == "categorical", field_code]
  filePath = sprintf(paste0(OUTPUT_PATH, "/%s/%s_%s_auprc.rds"), toupper(gene), gene, catePhenotypes)
  
  if (length(filePath) > 0 & all(file.exists(filePath))) {
    # Load AUPRCs
    results = lapply(filePath, readRDS)
    
    # Collect auprcs
    phenotypes = bPhenotypes[field_type == "categorical", as.character(field_code)]
    auprcs = lapply(1:length(results), function(index) {
      res = results[[index]]
      vals = res$auprcs
      vals$rn = sapply(vals$rn, function(val) names(VARIANT_PREDICTORS[VARIANT_PREDICTORS == val]), USE.NAMES = F)
      vals$code = phenotypes[[index]]
      vals$phenotype = bPhenotypes[field_id == phenotypes[[index]], field_description]
      return(vals[, .(predictor = rn, phenotype, code, auprc = auprc_avg, auprc_ci)])
    })
    auprcs = rbindlist(auprcs)
    
    # Collect p values
    prcPvals = lapply(results, "[[", "auprcs_pval")
    
    # Helper function: attach AUPRC to predictor name
    attachAUPRC = function(name, fieldId) {
      formName = unlist(VARIANT_PREDICTORS)
      formName = names(formName)[formName == name]
      val = auprcs[predictor == formName & code == fieldId, auprc]
      return(sprintf("%s (%.3f)", formName, val))
    }
    
    # Attach information
    prcPvals = mapply(function(tab, code) {
      tab$code = code
      tab$type = "PRC"
      tab$pred1 = sapply(tab$pred1, attachAUPRC, code)
      tab$pred2 = sapply(tab$pred2, attachAUPRC, code)
      return(tab)
    }, prcPvals, catePhenotypes, SIMPLIFY = F)
    prcPvals = rbindlist(prcPvals)
  }
  
  # Merge p values
  pVals = rbind(pccPvals, prcPvals)
  if (is.null(pVals)) return(NULL)
  
  # Compute q values
  phenotypes_ = unique(pVals$code)
  for (i in phenotypes_) {
    pVals[pVals$code == i, 'q_val'] = qvalue(pVals[pVals$code == i]$p_val)$qvalues
  }
  
  # Compute q-values (--previous method--)
  # q-values were being calculated based on the distribution
  # of all p-values for a given gene (rather than each gene-trait combination)
  #pVals$q_val = qvalue(pVals$p_val)$qvalues
  
  # Plot p-value heatmaps
  plots = lapply(unique(pVals$code), function(phenotype) {
    # Format labels
    plotTable = pVals[code == phenotype]
    plotTable$pred1 = unlist(plotTable$pred1)
    plotTable$pred2 = unlist(plotTable$pred2)
    type = unique(plotTable$type)
    
    # Plot matrix as heatmap
    plot = plotTable %>%
      dplyr::mutate(pred1 = forcats::fct_reorder(pred1, extractAUPRC(pred1)),
                    pred2 = forcats::fct_reorder(pred2, extractAUPRC(pred2))) %>%
      ggplot() +
      geom_tile(aes(x = pred2, y = pred1, fill = `p_val`), size = 0.5, color = "#d8d8d8") +
      scale_fill_gradientn(name = "P Val", limits = c(0, 0.05), labels = c("0", "\u2265 0.05"),
                           breaks = c(0, 0.05), na.value = "white",
                           colors = c("#034c9d", "white")) +
      theme_minimal(base_size = 16) +
      ggtitle(bPhenotypes[field_code == phenotype, field_description]) +
      labs(x = sprintf("Predictor Name (Mean %s)", type),
           y = sprintf("Predictor Name (Mean %s)", type)) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            legend.position = "bottom", plot.title = element_markdown(),
            legend.title = element_text(vjust = 0.75))
    
    return(plot)
  })
  
  # Arrange and save plots
  plots = arrangeGrob(grobs = plots, ncol = min(length(plots), 2))
  outputFile = sprintf(paste0(OUTPUT_PATH, "/%s/%s_pval_heatmaps.png"), toupper(gene), gene)
  ggsave(outputFile, plots, height = ceiling(length(plots) / 2) * 9,
         width = 9 * min(2, length(plots)), units = "in")
  
  # Plot q-value heatmaps
  plots = lapply(unique(pVals$code), function(phenotype) {
    # Format labels
    plotTable = pVals[code == phenotype]
    plotTable$pred1 = unlist(plotTable$pred1)
    plotTable$pred2 = unlist(plotTable$pred2)
    type = unique(plotTable$type)
    
    # Plot matrix as heatmap
    plot = plotTable %>%
      dplyr::mutate(pred1 = forcats::fct_reorder(pred1, extractAUPRC(pred1)),
                    pred2 = forcats::fct_reorder(pred2, extractAUPRC(pred2))) %>%
      ggplot() +
      geom_tile(aes(x = pred2, y = pred1, fill = `q_val`), size = 0.5, color = "#d8d8d8") +
      scale_fill_gradientn(name = "Q Val", limits = c(0, 0.1), labels = c("0", "\u2265 0.1"),
                           breaks = c(0, 0.1), na.value = "white",
                           colors = c("#034c9d", "white")) +
      theme_minimal(base_size = 16) +
      ggtitle(bPhenotypes[field_code == phenotype, field_description]) +
      labs(x = sprintf("Predictor Name (Mean %s)", type),
           y = sprintf("Predictor Name (Mean %s)", type)) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            legend.position = "bottom", plot.title = element_markdown(),
            legend.title = element_text(vjust = 0.75))
    
    return(plot)
  })
  
  # Arrange and save plots
  plots = arrangeGrob(grobs = plots, ncol = min(length(plots), 2))
  outputFile = sprintf(paste0(OUTPUT_PATH, "/%s/%s_qval_heatmaps.png"), toupper(gene), gene)
  ggsave(outputFile, plots, height = ceiling(length(plots) / 2) * 9,
         width = 9 * min(2, length(plots)), units = "in")
  
  # Add information and return
  pVals$gene = gene
  
  cat(sprintf("\rProcess gene %s (%d/%d)... Completed", gene, index, length(genes)), rep(" ", 15), "\n")
  flush.console()
  
  return(pVals)
})
pVals = rbindlist(pVals)
pVals = unique(pVals)

# Save p-value
fwrite(pVals, "pvals.csv")

# Helper function: extract AUPRC to predictor label
extractNumber = function(labels) {
  return(as.numeric(str_extract_all(labels, "\\d+\\.\\d+|-\\d+\\.\\d+")))
}

# Extract score
pVals$pred1_score = extractNumber(pVals$pred1)
pVals[type == "PCC", pred1_score := pred1_score^2]
pVals$pred2_score = extractNumber(pVals$pred2)
pVals[type == "PCC", pred2_score := pred2_score^2]
pVals$pred1_name = str_split(pVals$pred1, " ", simplify = T)[, 1]
pVals$pred2_name = str_split(pVals$pred2, " ", simplify = T)[, 1]

# Restrict to gene-trait combinations where every predictor was able to make prediction
#subsetGeneCombs = pVals[, .N, by = c("gene", "code")][N == 441] # 21 * 21 = 441 pairwise predictor comparisons
#pVals = merge(pVals, subsetGeneCombs[, .(gene, code)])

# Count the number of variant predictors that are indistinguishable from 
# the best performing predictor
predictors = apply(pVals[, .SD[which.max(pred1_score)], by = c("gene", "code")], 1, function(row) {
  # Select relevant variant predictors
  preds = pVals[gene == row[["gene"]] & pred1 == row[["pred1"]]]
  
  # Select indistinguishable preds
  preds = preds[q_val >= 0.10, pred2]
  
  if (length(preds) < 1) return(NULL)
  
  return(data.table(gene = row[["gene"]],
                    indistingushable_pred = str_split(preds, fixed(" ("), simplify = T)[, 1]))
})
predictors = rbindlist(predictors)

# Summarize predictors
predictors = predictors[, .(num_indistinguishable = .N), by = "indistingushable_pred"]

# Compute Wilcoxon U test p values
sigTest = pVals[, .N, by = c("pred1_name", "pred2_name")]
sigs = apply(sigTest, 1, function(row) {
  pred1Name = row[["pred1_name"]]
  pred2Name = row[["pred2_name"]]
  
  # Do not compare two same predictors
  if (pred1Name == pred2Name) return(NA)
  
  # Do not compare if predictor 2 did not perform as well as predictor 1
  pred1Perf = predictors[indistingushable_pred == pred1Name, num_indistinguishable]
  pred2Perf = predictors[indistingushable_pred == pred2Name, num_indistinguishable]
  if (pred2Perf <= pred1Perf) return(1)
  
  # Filtered table
  subPVals = pVals[pred1_name == pred1Name & pred2_name == pred2Name]
  
  sig = wilcox.test(subPVals$pred1_score, subPVals$pred2_score, alternative = "two.sided", paired = T, exact = F)
  return(sig$p.value)
})
sigTest$p_val = sigs
sigTest$q_val = qvalue(sigs)$qvalues

# Plot comparison
plotTable = merge(sigTest, predictors, by.x = "pred1_name", by.y = "indistingushable_pred")
plotTable = merge(plotTable, predictors, by.x = "pred2_name", by.y = "indistingushable_pred")
plot = plotTable %>%
  dplyr::mutate(pred1 = sprintf("%s (%d)", pred1_name, num_indistinguishable.x),
                pred2 = sprintf("%s (%d)", pred2_name, num_indistinguishable.y),
                q_val = ifelse(q_val > 0.1, 0.1, q_val)) %>%
  dplyr::mutate(pred1 = forcats::fct_reorder(pred1, num_indistinguishable.x),
                pred2 = forcats::fct_reorder(pred2, num_indistinguishable.y)) %>%
  ggplot() +
  geom_tile(aes(x = pred1, y = pred2, fill = q_val), linewidth = 0.5, color = "#d8d8d8") +
  scale_fill_gradientn(name = "FDR", limits = c(0, 0.1), labels = c("0", "\u2265 0.1"),
                       breaks = c(0, 0.1), na.value = "lightgrey",
                       colors = c("#034c9d", "white")) +
  theme_minimal(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = "bottom", plot.title = element_markdown(),
        legend.title = element_text(vjust = 0.75))
ggsave("method_comparsion.png", plot, height = 7, width = 7, units = "in")

# Plot absolute indistingusiable
plotTable = predictors[, .(method = indistingushable_pred, num_indistinguishable)]
plot = ggplot(plotTable, aes(x = reorder(method, num_indistinguishable), 
                             y = num_indistinguishable)) +
  geom_bar(stat="identity") +
  geom_text(aes(label = num_indistinguishable), hjust = 1.5, colour = "white") +
  labs(x = "Variant Effect Predictor",
       y = "Number of phenotype-gene combinations\nwhere the predictor is indistinguishable from the best performing predictor") +
  ggtitle(sprintf("Number of phenotype-gene combinations considered: %d", 
                  nrow(pVals[, .N, by = c("gene", "code")]))) +
  coord_flip() + theme_minimal(base_size = 14)
outputFile = sprintf("predictors_benchmark.png")
ggsave(outputFile, plot, width = 10, height = 8, units = "in")
