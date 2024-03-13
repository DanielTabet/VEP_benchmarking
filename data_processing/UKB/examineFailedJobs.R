library(data.table)
library(stringr)

# Load global variables from the command line
args = commandArgs(trailingOnly = T)
TOTAL_JOB_PATH = args[1]
FAILED_JOB_PATH = args[2]
VCF_FILE_PATH = args[3]
OUTPUT_VCF_PATH = args[4]

# TOTAL_JOB_PATH = "submitted_jobs_parse_vep_chr10-19.log"
# FAILED_JOB_PATH = "failed_jobs_parse_vep_chr1-19.txt"
# VCF_FILE_PATH = "/Bulk/Exome sequences/Population level exome OQFE variants, pVCF format - interim 450k release"
# OUTPUT_VCF_PATH = "vcf_files_failed_parse_vep_chr10-19.txt"

# Load total jobs
totalJobs = readLines(TOTAL_JOB_PATH)
totalJobs = str_split(totalJobs, fixed("; "), simplify = T)
totalJobs[, 1] = str_split(totalJobs[, 1], fixed(": "), simplify = T)[, 2]
totalJobs[, 2] = str_split(totalJobs[, 2], fixed(": "), simplify = T)[, 2]
totalJobs = as.data.table(totalJobs)
colnames(totalJobs) = c("vcf_file", "job_id")

# Load failed jobs
failedJobs = readLines(FAILED_JOB_PATH)

# Filter total jobs to find failed jobs with associated VCF file
totalJobs[, failed := F]
totalJobs[job_id %in% failedJobs, failed := T]

# Save failed input VCF files for resubmission
writeLines(totalJobs[failed == T, vcf_file], OUTPUT_VCF_PATH)