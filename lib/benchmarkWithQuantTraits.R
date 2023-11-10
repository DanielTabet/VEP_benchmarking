library(data.table)
library(ggplot2)
library(stringr)
library(ggtext)
library(forcats)
library(Rmisc)
library(gridExtra)
library(fs)
library(dplyr)
library(grid)

# Set arguments if not passed in already as config
GENE = tolower(config$gene)
VARIANT_PREDICTORS = config$variant_predictors
PLOT_INDIVIDUAL_CORRELATIONS = config$plot_individual_correlations
BOOTSTRAP_N = config$bootstrap_iterations
VARIANTS_CUTOFF = config$occurance_cutoff
OUTPUT_PATH = config$output_path

###
# Process variants with quantitative phenotypes
# Stats: PCC
###
# Load variants with quantitative phenotypes
filePath = sprintf(paste0(OUTPUT_PATH, "/%s/%s_quant-variants_measures.csv"), toupper(GENE), GENE)
if (!file.exists(filePath)) stop("no quant variants")
quantVariants = fread(filePath, header = T)

# Load variant prediction (VP) scores
varPredictors = quantVariants[, c("variant_name", unlist(VARIANT_PREDICTORS, use.names = F)), with = F]
varPredictors = unique(varPredictors)

# Remove predictors that didn't make predictions for at least the number of variants set by the VARIANTS_CUTOFF
colsKept = apply(varPredictors[, -c("variant_name"), with = F], 2, function(col) {
  length(na.omit(col)) >= VARIANTS_CUTOFF
})
if (sum(colsKept) < 1) stop("no variant effect predictors made enough predictions")
varPredictors = varPredictors[, -names(colsKept[!colsKept]), with = F]
cols = names(colsKept[colsKept])
VARIANT_PREDICTORS = VARIANT_PREDICTORS[VARIANT_PREDICTORS %in% cols]

# Helper function: calculate PCC and plot (optional)
calculateCorrAndPlot = function(code, type, variants, returnPlot = F) {
  description = bPhenotypes[field_id == code, field_description]
  results = lapply(VARIANT_PREDICTORS, function(predictor) {
    plotTable = variants[, c(predictor, "p_score"), with = F]
    plotTable = plotTable[complete.cases(plotTable)]
    setnames(plotTable, predictor, "vp_score")
    varPredLabel = names(VARIANT_PREDICTORS[VARIANT_PREDICTORS == predictor])
    corrs = c("PCC" = cor(plotTable$vp_score, plotTable$p_score, method = "pearson"))
    
    if (returnPlot) {
      plot = ggplot(plotTable, aes(x = p_score, y = vp_score)) +
        geom_point() +
        geom_smooth(method='lm', formula= y~x, se = F) +
        labs(x = description, y = paste(varPredLabel, "score")) +
        theme_minimal(base_size = 16) +
        ggtitle(sprintf("%s - %s", toupper(GENE), description),
                sprintf("PCC = %.2f, N = %d", corrs[[1]], nrow(plotTable)))
      
      return(list("name" = varPredLabel, "corrs" = corrs, "plot" = plot))
    }
    
    return(list("name" = varPredLabel, "corrs" = corrs))
  })
  
  # Collect correlations
  correlations = lapply(results, function (result) {
    return(data.table(phenotype = description, name = result$name, type = type,
                      PCC = result$corrs[["PCC"]]))
  })
  correlations = rbindlist(correlations)
  
  # Merge plots
  if (PLOT_INDIVIDUAL_CORRELATIONS && returnPlot) {
    plots = lapply(results, "[[", "plot")
    plots = arrangeGrob(grobs = plots, ncol = 4)
    
    # Save plots
    plotPath = sprintf(paste0(OUTPUT_PATH, "/%s/%s_%s.png"), toupper(GENE), GENE, path_sanitize(description))
    ggsave(plotPath, plots, dpi = 300,
           units = "in", width = 20, height = 5 * ceiling(length(plots) / 4))
  }
  
  return(correlations)
}

# Helper function: calculate std err
calcStdErr = function(x) sd(x)/sqrt(length(x))

# Helper function: format numbers with scientific notation
changeSciNot = function(n) {
  output = format(n, scientific = TRUE) #Transforms the number into scientific notation even if small
  output = sub("e", "*10^", output) #Replace e with 10^
  output = sub("\\+0?", "", output) #Remove + symbol and leading zeros on expoent, if > 1
  output = sub("-0?", "-", output) #Leaves - symbol but removes leading zeros on expoent, if < 1
  output
}

# Helper function: benchmark phenotypes
benchmarkPhenotype = function(phenotype) {
  # Extract phenotype code, type and description
  code = str_trim(as.character(phenotype[["field_id"]]))
  type = phenotype[["field_type"]]
  description = phenotype[["field_description"]]
  
  # Select variants based on phenotype type
  # Use deep copy to make sure the original data table is not changed
  variants = copy(quantVariants)
  if (!code %in% colnames(variants))
    stop(sprintf("Missing measurement for phenotype %s", code))
  setnames(variants, code, "p_score")
  
  # Use bootstrap resampling to determine variants correlations and CI
  numCount = 0
  correlations = data.table()
  
  while (numCount < BOOTSTRAP_N) {
    # Sample variants
    sampleVariants = variants[sample(1:nrow(variants), nrow(variants), replace = T),
                              .(variant_name, p_score)]
    
    # Calculate mean phenotype score of the samples
    sampleVariants = sampleVariants[, .(p_score = mean(p_score, na.rm = T)), by = "variant_name"]
    sampleVariants = sampleVariants[complete.cases(sampleVariants)]
    
    # Attach VP scores to phenotype
    sampleVariants = merge(sampleVariants, varPredictors, by = "variant_name")
    
    # Calculate and return correlations
    corr = calculateCorrAndPlot(code, type, sampleVariants)
    
    # Process sample result
    numCount = numCount + 1
    cat("\rBenchmarking phenotype", description, "with bootstrapping, iteration", numCount)
    flush.console()
    if (nrow(correlations) < 1) {
      # If correlations are empty (just created), assign the results of the 
      # first run
      correlations = corr
    } else {
      # Otherwise
      # Add sample results to list
      correlations = merge(correlations, corr, by = c("phenotype", "name", "type"))
      correlations$PCC = mapply(c, correlations$PCC.x, correlations$PCC.y, SIMPLIFY = F)
      correlations = correlations[, .(phenotype, name, type, PCC)]
    }
  }
  
  # Print complete message
  cat("\rBenchmarking phenotype", description, "with bootstrapping... Completed",
      rep(" ", 15), "\n")
  flush.console()
  
  # Plot bootstrap-averaged and empirical correlations
  results = lapply(1:nrow(correlations), function(index) {
    # Get basic info
    phenotype = correlations[index, phenotype]
    predictor = correlations[index, name]
    
    # Check for NA
    plotTable = correlations[index, .(pcc = na.omit(unlist(PCC)))]
    if (nrow(plotTable) < 1) {
      ret = correlations[name == predictor, .(phenotype, name, type, pcc = PCC,
                                              avg_pcc = NA, avg_pcc_ci = NA)]
      return(list(plots = NA, correlations = ret))
    }
    
    # Plot average PCC distribution
    pccCi = quantile(plotTable$pcc, c(0.05, 0.95))
    distribution = ggplot(plotTable, aes(x = pcc)) +
      geom_histogram(position="identity", 
                     alpha = 0.5, bins = 50, fill = "#f46d43") +
      ggtitle(sprintf("%s - %s", phenotype, predictor),
              sprintf("Final Average PCC = %.3f, Stdev = %.3f, 5-95%% Range = [%.3f,%.3f]",
                      mean(plotTable$pcc), sd(plotTable$pcc), 
                      pccCi[["5%"]], pccCi[["95%"]])) +
      labs(x = "Correlation", y = "Bootstrap Sample Count") +
      theme_minimal(base_size = 16)
    
    # Return bootstrap and empirical correlations and 95% CI
    ret = correlations[name == predictor, .(phenotype, name, type, pcc = PCC,
                                            avg_pcc = mean(plotTable$pcc),
                                            avg_pcc_ci = list(pccCi[c("5%", "95%")]))]
    return(list(plots = distribution, correlations = ret))
  })
  
  # Collect correlations
  correlations = rbindlist(lapply(results, "[[", 2))
  correlations$code = code
  
  # Merge plots and save
  plots = results[!sapply(results, anyNA, recursive = F)]
  plots = arrangeGrob(grobs = lapply(plots, "[[", 1), ncol = 4,
                      top = textGrob(paste("Bootstrapping N =", BOOTSTRAP_N),
                                     gp = gpar(fontsize = 20)))
  ggsave(sprintf(paste0(OUTPUT_PATH, "/%s/%s_%s_bootstrap.png"), toupper(GENE), GENE, path_sanitize(description)),
         plots, dpi = 300, units = "in",
         width = floor(length(results) / 4) * 7,
         height = floor(length(results) / 4) * 5)
  
  return(correlations)
}

# Benchmark quantitative phenotypes
if (nrow(bPhenotypes[field_type == "quantitative"]) < 1) {
  cat("No quantitative phenotypes")
} else {
  correlations = apply(bPhenotypes[field_type == "quantitative"],
                       1, benchmarkPhenotype)
  correlations = rbindlist(correlations)
  correlations$index = 1:nrow(correlations)
  
  # Helper function: format CI
  formatCI = function(range) {
    if (any(is.na(range))) return("NA")
    
    return(paste(formatC(round(range, 2), format='f', digits=2), collapse = ",<br>"))
  }
  
  # Format CI
  plotTable = correlations[, .(phenotype, predictor = name, type, index, avg_score = avg_pcc,
                               avg_score_ci = lapply(avg_pcc_ci, formatCI))]
  
  # Order predictors by its average performance
  predictorPerf = plotTable[, .(avg_score = mean(abs(avg_score))), by = "predictor"]
  predictorPerf[order(avg_score, decreasing = T), predictor_index := 1:nrow(predictorPerf)]
  plotTable = merge(plotTable, predictorPerf[, .(predictor, predictor_index)], by = "predictor")
  
  # Plot PCC heatmap
  plot = plotTable %>%
    mutate(phenotype = sprintf("%s<br>*(%s)*", phenotype, type),
           phenotype = fct_reorder(phenotype, -index),
           name = fct_reorder(predictor, -predictor_index)) %>%
    ggplot(aes(x = name, y = phenotype, fill = avg_score)) +
    geom_tile() +
    geom_richtext(aes(label = sprintf("%.2f<br>[%s]", avg_score, avg_score_ci)),
                  size=5, fill = NA, label.color = NA) +
    labs(x = "Variant Effect Predictors<br><span style='font-size:12pt'>(ordered right-to-left descendingly by the mean of absolute scores)</span>",
         y = "Phenotype",
         fill = "Mean PCC") +
    scale_fill_gradientn(limits = c(-1, 1), breaks = c(-1, 0, 1),
                         colors = c("#ef8a62", "#ffffff", "#67a9cf"), na.value = "#d8d8d8") +
    ggtitle(paste(toupper(GENE), "with bootstrapping, N =", BOOTSTRAP_N),
            "5-95% range in bracket") +
    theme_minimal(base_size = 18) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 20),
          axis.text.y = element_markdown(size = 20), legend.position = "bottom",
          axis.title.x = element_markdown(),
          legend.title.align = 1, legend.title = element_markdown(size = 16, vjust = 0.9))
  
  ggsave(sprintf(paste0(OUTPUT_PATH, "/%s/%s_quant_bootstrap.png"), toupper(GENE), GENE), plot, dpi = 300, 
         height = 4.5 + length(unique(plotTable$phenotype)), 
         width = max(10, floor(length(VARIANT_PREDICTORS) * 1.1)), units = "in")
  
  # Format and save scores
  setnames(correlations, c("phenotype", "code"), c("field_description", "field_id"))
  outputFile = sprintf(paste0(OUTPUT_PATH, "/%s/%s_correlations.csv"), toupper(GENE), GENE)
  fwrite(correlations, outputFile)
}