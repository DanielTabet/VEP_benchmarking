library(data.table)
library(stringr)

# Set arguments
GENE = tolower(config$gene)
ENSEMBL_ID = toupper(config$ensembl_id)
#CAN_TRANSCRIPT_ID = toupper(config$canonical_transcript_id)
INPUT_VARIANTs_DIR = config$input_var_dir
ALL_VARIANTs_PATH = config$all_variants_path
UNIQUE_VARIANTs_PATH = config$unique_variants_path
WITHDRAW_EIDs_PATH = config$withdraw_eids_path
VARIANT_BLOCKs_PATH = config$variant_blocks_path
OUTPUT_PATH = config$output_path

# SELECT VARIANTS FOR GENE OF INTEREST
# 1) Get genomic coordinates
ensemblDB = fread("common/ensemblDB.txt")
coordinate = data.table(ensemblDB)
coordinate = coordinate[ensemblID==ENSEMBL_ID]

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

# SELECT VARIANTS
# 1) Allele frequency -- gnomAD MAF < 0.1% (0.001) | UKB < 0.1% (0.001)
mergedVariants[is.na(gnomAD_AF), gnomAD_AF := 0]
mergedVariants = mergedVariants[gnomAD_AF < 0.001 & AF < 0.001]

# 2) Missense variants
mergedVariants = mergedVariants[str_detect(Consequence, "missense_variant")]

# 3) Clean up variants
mergedVariants = mergedVariants[str_detect(Protein_position, "-", negate=TRUE)]
mergedVariants = mergedVariants[str_detect(CDS_position, "-", negate=TRUE)]

# 4) Select columns
mergedVariants =  mergedVariants[, c('variant_name', 'eid', 'chr', 'pos', 'ref', 'alt', 'Gene',
                                     'Consequence', 'cDNA_position', 'CDS_position', 'Protein_position',
                                     'Amino_acids', 'Codons', 'gnomAD_AF', 'AF')]

# LOAD VARIANT EFFECT PREDICTIONS
# Column list
score_cols = c('UniProtID',	'aapos', 'aaref',	'aaalt',	'Ensembl_geneid',	'Ensembl_transcriptid',
               'variant_name',	'REVEL_score',	'VARITY_R_score',	'Polyphen2_HVAR_score',
               'PROVEAN_score',	'SIFT_score',	'FATHMM_score',	'MPC_score',	'LRT_score',	
               'PrimateAI_score',	'CADD_raw',	'DANN_score', 'Eigen-raw_coding',	'GenoCanyon_score',
               'M-CAP_score',	'MetaLR_score',	'MetaSVM_score',	'MVP_score',	'MutationTaster_score',
               'SiPhy_29way_logOdds',	'alphaMissense_score',	'EVE_score',	'MutPred2_score',
               'ESM1b_score',	'ESM1v_score')

# VEPs pre-processed dbNSFPv4 file with add-ons
scores = fread('common/VEP_scores_vFinal.csv', select=score_cols)

# Merge predictor scores for gene with UKB variants
mergedVariants = merge(mergedVariants, scores[Ensembl_geneid == ENSEMBL_ID],
                       by.x = 'variant_name', by.y = 'variant_name')

# Flip some predictor scores
# PROVEAN, SIFT, FATHMM, LRT, ESM1b and ESM1v
mergedVariants[, PROVEAN_flipped := -PROVEAN_score]
mergedVariants[, SIFT_flipped := -SIFT_score]
mergedVariants[, FATHMM_flipped := -FATHMM_score]
mergedVariants[, LRT_flipped := -LRT_score]
mergedVariants[, ESM1b_flipped := -ESM1b_score]
mergedVariants[, ESM1v_flipped := -ESM1v_score]

# Add CADD2 scores
# These were only calculated for the UKB set and are not in the pre-processed file
CADD2 = fread("common/CADD2/FritzDaniel_LR-Ce-01-wmo-i013.tsv",
              select=c('oAA', 'nAA', 'protPos', 'GeneID', 'RawScore'))

# Merge
geneID = unique(mergedVariants$Ensembl_geneid)
mergedVariants = merge(mergedVariants, CADD2[GeneID = geneID],
                       by.x = c("Ensembl_geneid", "aapos", "aaref", "aaalt"),
                       by.y = c("GeneID", "protPos", "oAA", "nAA"), all.x=T)

# Select unique entry ID (eid)
eids = sort(unique(mergedVariants$eid))

# Remove withdrawn participants
eids = eids[eids > 0] # Negative EID are assigned to withdrawn participants
withdraws = readLines(WITHDRAW_EIDs_PATH)
withdraws = as.numeric(withdraws)
eids = eids[!eids %in% withdraws]

# Save unique entry ID and the merged variants
outputPrefix = sprintf(paste0(OUTPUT_PATH, "/%s/%s"), toupper(GENE), GENE)
writeLines(as.character(eids), paste0(outputPrefix, "unique_eids_filtered.txt"))
fwrite(mergedVariants, paste0(outputPrefix, "variants_filtered.csv"))
