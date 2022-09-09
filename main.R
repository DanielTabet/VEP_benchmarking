source("lib/packageCheck.R")

# Check argument
args = commandArgs(trailingOnly = T)
if (length(args) != 2) {
  cat("Usage:\n")
  cat("  Rscript main.R <configuration-file> <log-dir>")
  quit()
}

# Load configuration file
source("lib/loadConfig.R")
config = loadConfig(args[1])

# Start logging
if (is.na(args[2])) stop("Log dir path not specified")
library(logr)
logFile = log_open(paste(args[2], "main.log", sep = "/"))

# If no output folder, create one
if (!dir.exists("output")) dir.create("output")

# Process gene list
library(data.table)
geneList = fread(config$gene_list)

# Validate that gene inputs are available
hasInput = dir.exists(paste("input", geneList$gene_symbol, sep = "/"))
if (!all(hasInput)) {
  warning("Missing gene input for genes:",
       paste0(geneList$gene_symbol[!hasInput], collapse = ", "))
}

# Start workers
invisible(apply(geneList, 1, function (row) {
  gene = row[["gene_symbol"]]
  ensemblId = row[["gene_id"]]
  transcriptId = row[["canonical_transcript"]]
  index = which(geneList$gene_symbol == gene)
  log_print(paste("==> Process gene", gene,
                  sprintf("(index = %d)", index)))
  logPath = sprintf("%s.log", gene)
  system2(config$rscript_path, args = c("worker.R", args[1], gene, ensemblId,
                                        transcriptId,
                                        paste(args[2], logPath, sep = "/")))
  log_print("")
}))

# Collect benchmark results
log_print("Collecting benchmark results")
source("lib/collectBenchmarkResults.R")

# Compute misc stats
log_print("Computing misc. stats")
source("lib/miscStats.R")

log_print("Complete!")
log_close()