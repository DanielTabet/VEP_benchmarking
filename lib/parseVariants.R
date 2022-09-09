library(data.table)
library(stringr)
library(EnsDb.Hsapiens.v86)
Sys.setenv(PATH=paste(Sys.getenv("PATH"), "/opt/homebrew/bin/", sep=":"))

# Set arguments if not passed in already as config
GENE = tolower(config$gene)
ENSEMBL_ID = toupper(config$ensembl_id)
CAN_TRANSCRIPT_ID = toupper(config$canonical_transcript_id)
INPUT_VARIANTs_DIR = config$input_var_dir
ALL_VARIANTs_PATH = config$all_variants_path
UNIQUE_VARIANTs_PATH = config$unique_variants_path
WITHDRAW_EIDs_PATH = config$withdraw_eids_path
VARIANT_BLOCKs_PATH = config$variant_blocks_path

# Identify variants belong to the gene of interest
# 1) Get gene's genome coordinate
edb = EnsDb.Hsapiens.v86
gene = genes(edb, filter = GeneIdFilter(ENSEMBL_ID))
coordinate = data.table("chr" = as.character(seqnames(gene)),
                        "start" = start(gene), "end" = end(gene))
# 2) Select variant files within range
varBlocks = fread(VARIANT_BLOCKs_PATH, drop = 1, col.names = c("chr", "block", "start", "end"))
blocks = varBlocks[chr == coordinate$chr]
blocks = blocks[start <= coordinate$start & end >= coordinate$end]
if (!nrow(blocks)) {
  # The block spans over multiple blocks
  # handle this case separately
  blocks = varBlocks[chr == coordinate$chr]
  blocks = rbind(blocks[start <= coordinate$start][which.max(start)],
                 blocks[end >= coordinate$end][which.min(start)])
}
# If still no blocks, stop
if (!nrow(blocks)) stop("no blocks found. Please check the genome coordinate")
# 3) Load variants within the selected blocks
mergedVariants = apply(blocks, 1, function(row) {
  chrom = row[["chr"]]
  block = row[["block"]]
  
  # Load variants
  fileName = sprintf(ALL_VARIANTs_PATH, chrom, block)
  allVariants = fread(paste(INPUT_VARIANTs_DIR, fileName, sep = "/"), drop = 1)
  fileName = sprintf(UNIQUE_VARIANTs_PATH, chrom, block)
  uniqueVariants = fread(paste(INPUT_VARIANTs_DIR, fileName, sep = "/"))
  
  # Subset all variants
  mergedVariants = merge(allVariants, uniqueVariants[Gene == ENSEMBL_ID],
                         by.x = "variant_name", by.y = "#Uploaded_variation")
  
  return(mergedVariants)
})
mergedVariants = rbindlist(mergedVariants, use.names = T, fill = T)

# Select unique variants:
# 1) gnomAD MAF < 0.1% (0.001)
# 2) Missense
mergedVariants = mergedVariants[is.na(gnomAD_AF), gnomAD_AF := 0]
mergedVariants = mergedVariants[gnomAD_AF < 0.001 &
                                  str_detect(Consequence, "missense_variant")]

# Split scores
mergedVariants[, Ensembl_transcriptid := str_split(Ensembl_transcriptid, fixed(","))]
mergedVariants[, Uniprot_acc := str_split(Uniprot_acc, fixed(","))]
mergedVariants[, PROVEAN_score := str_split(PROVEAN_score, fixed(","))]
mergedVariants[, SIFT_score := str_split(SIFT_score, fixed(","))]
mergedVariants[, FATHMM_score := str_split(FATHMM_score, fixed(","))]
mergedVariants[, Polyphen2_HVAR_score := str_split(Polyphen2_HVAR_score, fixed(","))]
mergedVariants[, MPC_score := str_split(MPC_score, fixed(","))]
mergedVariants[, MVP_score := str_split(MVP_score, fixed(","))]
mergedVariants[, MutationTaster_score := str_split(MutationTaster_score, fixed(","))]
mergedVariants = suppressWarnings(lapply(1:nrow(mergedVariants), function(index) {
  row = mergedVariants[index,]
  transcripts = unlist(row$Ensembl_transcriptid)
  transIndex = which(CAN_TRANSCRIPT_ID == transcripts)
  # If nothing found, use the first transcript
  if (length(transIndex) < 1) {
    warning("mismatch transcript. Use the first transcript.")
    transIndex = 1
  }
  row$Uniprot_acc = unlist(row$Uniprot_acc)[transIndex]
  row$PROVEAN_score = as.numeric(unlist(row$PROVEAN_score)[transIndex])
  row$SIFT_score = as.numeric(unlist(row$SIFT_score)[transIndex])
  row$FATHMM_score = as.numeric(unlist(row$FATHMM_score)[transIndex])
  row$Polyphen2_HVAR_score = as.numeric(unlist(row$Polyphen2_HVAR_score)[transIndex])
  row$MPC_score = as.numeric(unlist(row$MPC_score)[transIndex])
  row$MVP_score = as.numeric(unlist(row$MVP_score)[transIndex])
  row$MutationTaster_score = as.numeric(unlist(row$MutationTaster_score)[transIndex])
  return(row)
}))
mergedVariants = rbindlist(mergedVariants)

# Flip some predictor scores
# provean, sift, fathmm, LRT
mergedVariants[, PROVEAN_flipped := -PROVEAN_score]
mergedVariants[, SIFT_flipped := -SIFT_score]
mergedVariants[, FATHMM_flipped := -FATHMM_score]
mergedVariants[, LRT_flipped := -LRT_score]

# Merge VARITY
mergedVariants$VARITY_ER = NULL
mergedVariants$VARITY_R = NULL
mergedVariants$Uniprot_acc = str_split(mergedVariants$Uniprot_acc, fixed("-"), simplify = T)[,1]
mergedVariants$Protein_position = as.numeric(mergedVariants$Protein_position)
varity = fread("common/varity_all_predictions.txt")
mergedVariants = merge(mergedVariants,
                       varity[, .(p_vid, nt_alt, aa_pos, aa_ref, aa_alt, VARITY_R, VARITY_ER, VARITY_R_LOO, VARITY_ER_LOO)],
                       by.x = c("Uniprot_acc", "alt", "Protein_position", "aaref", "aaalt"),
                       by.y = c("p_vid", "nt_alt", "aa_pos", "aa_ref", "aa_alt"),
                       all.x = T)

# Merge EVE
eve = fread("common/eve_variants.csv")
mergedVariants = merge(mergedVariants,
                       eve[, .(uniprot_id, wt_aa, position, mt_aa, EVE_score = EVE_scores_ASM)], 
                       by.x = c("Uniprot_acc", "Protein_position", "aaref", "aaalt"),
                       by.y = c("uniprot_id", "position", "wt_aa", "mt_aa"),
                       all.x = T)

# Select unique entry ID (eid)
eids = sort(unique(mergedVariants$eid))

# Remove withdrawn participants
eids = eids[eids > 0] # Negative EID are assigned to withdrawn participants
withdraws = readLines(WITHDRAW_EIDs_PATH)
withdraws = as.numeric(withdraws)
eids = eids[!eids %in% withdraws]

# Save unique entry ID and the merged variants
outputPrefix = sprintf("output/%s/%s_", toupper(GENE), GENE)
writeLines(as.character(eids), paste0(outputPrefix, "unique_eids_filtered.txt"))
fwrite(mergedVariants, paste0(outputPrefix, "variants_filtered.csv"))