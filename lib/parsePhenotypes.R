library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(tools)
library(reshape2)
library(odbc)
library(RMariaDB)
library(retry)

# Set arguments if not passed in already as config
GENE = tolower(config$gene)
DBCONNECTFILEPATH = config$db_connect_path

# Load database (DB) parameters
source("lib/loadConfig.R")
dbParams = loadConfig(DBCONNECTFILEPATH)

# Connect to DB
con <- dbConnect(
  drv = MariaDB(), 
  username = dbParams$username,
  password = dbParams$password, 
  host = dbParams$host,
  port = dbParams$port,
  dbname = dbParams$dbname
)

# Load participants
eids = readLines(sprintf("output/%s/%s_unique_eids_filtered.txt", toupper(GENE), GENE))
eids = paste("'", eids, "'", sep = "", collapse = ",")

# Query phenotypes
pids = bPhenotypes$field_id
pids = paste("'", pids, "'", sep = "", collapse = ",")
# Because in the database, pids are stored in the format of <phenotype>-<visit>-<measurement>,
# when matching phenotypes, we need to only match the first part (i.e. <phenotype>)
# to the phenotype IDs. This is why we need to first substring the pid in the SQL query below
query = sprintf("SELECT * FROM phenotypes WHERE eid IN (%s) AND SUBSTRING_INDEX(pid, '-', 1) IN (%s)",
                eids, pids)
res = dbGetQuery(con, query)
dbDisconnect(con) # Disconnect DB
phenotypes = as.data.table(res)
phenotypes$field_id = as.numeric(str_split(phenotypes$pid, "-", simplify = T)[,1])

# Update EID
eids = unique(phenotypes$eid)

# If field ID has underscore, it's a combinatorial phenotype
# that needs to be further filtered
if (any(str_detect(bPhenotypes$field_code, "_"))) {
  bPhenotypes[str_detect(field_code, "_"),
              secondary_field_code := str_split(field_code, "_", simplify = T)[, 2]]
  phenotypes = merge(phenotypes, bPhenotypes[, .(field_id, secondary_field_code)], by = "field_id")
  phenotypes = phenotypes[(!is.na(secondary_field_code) & str_starts(measurement, secondary_field_code))
                          | is.na(secondary_field_code)]
  phenotypes[!is.na(secondary_field_code), measurement := secondary_field_code]
}

# Remove negative measures to NA
phenotypes = phenotypes[measurement >= 0]

# Extract visit (first number after the dash)
phenotypes[, visit := as.numeric(str_sub(str_extract_all(pid, "-\\d", simplify = T), 2))]

# Only keep the most recent visit for each sample
phenotypes = phenotypes[, .SD[visit == max(visit)], by = c("eid", "field_id")]

# Helper function: select measurement
selectMeasurement = function(measurements) {
  # If numbers, get mean
  if (all(!is.na(suppressWarnings(as.numeric(measurements))))) {
    return(as.character(mean(as.numeric(measurements))))
  }
  
  # If not numbers, return unique
  return(unique(measurements))
}

# Make sure one phenotype measurement per eid
phenotypes = phenotypes[, .(measurement = selectMeasurement(measurement)), by = c("eid", "field_id")]
if (nrow(unique(phenotypes[, .(eid, field_id)])) != nrow(phenotypes)) {
  stop("some phenotypes are not unique for gene ", toupper(GENE))
}

# Reshape data table
phenotypes = acast(phenotypes, eid~field_id, value.var = "measurement")
phenotypes = as.data.table(phenotypes, keep.rownames = T)
setnames(phenotypes, "rn", "eid")

# Add missing EIDs
phenotypes = merge(data.table(eid = eids), phenotypes, all.x = T)

# Save filtered phenotypes
outputFile = sprintf("output/%s/%s_phenotypes_filtered.csv", 
                     toupper(GENE), GENE)
fwrite(phenotypes, outputFile)
