############################################################################
#This code provides sample scripts for processing phenotype data from All of Us (AoU)
#Used in Jupyter Notebooks in All of Us workbench (R kernel)
#This assumes that raw phenotype data has already been retrieved using the AoU-generated SQL and saved to your cloud storage bucket
#That process is outlined in AoU user resources. Phenotypes are defined using "concept sets."

#Raw phenotype data should contain columns: 
##'person_id' (unique subject id) 
##'unit_source_value' (units) 
##'value_as_number' (measurement value) 
##'measurement_datetime' (measurement data and time)
##'measurement_concept_id' (measurement code)
##'standard_concept_name' (measurement name)
##'source_concept_code'

############################################################################
#SETUP

##########
#libraries
library(tidyverse)

##########
#Define helper functions

#Retrieve raw data stored in Google Cloud Storage bucket
read_bq_export_from_workspace_bucket <- function(export_path) {
  col_types <- cols(standard_concept_name = col_character(), standard_concept_code = col_character(), standard_vocabulary = col_character(), observation_type_concept_name = col_character(), value_as_string = col_character(), value_as_concept_name = col_character(), qualifier_concept_name = col_character(), unit_concept_name = col_character(), visit_occurrence_concept_name = col_character(), observation_source_value = col_character(), source_concept_name = col_character(), source_concept_code = col_character(), source_vocabulary = col_character(), unit_source_value = col_character(), qualifier_source_value = col_character(), value_source_value = col_character())
  bind_rows(
    map(system2('gsutil', args = c('ls', export_path), stdout = TRUE, stderr = TRUE),
        function(csv) {
          message(str_glue('Loading {csv}.'))
          chunk <- read_csv(pipe(str_glue('gsutil cat {csv}')), col_types = col_types, show_col_types = FALSE)
          if (is.null(col_types)) {
            col_types <- spec(chunk)
          }
          chunk
        }))
}

############
#get path to Cloud Storage bucket associated with your workspace
your_bucket <- Sys.getenv('WORKSPACE_BUCKET')

#################################################################################
# CLEAN AND FORMAT QUANTITATIVE TRAITS

#Example: code for Hemoglobin A1c

#define path to raw data
A1c_path <- paste(your_bucket,'/data/your_trait_file_here.csv')

#define UK Biobank field id that corresponds to your trait (A1c in this example)
UKB_code <- 'code_30750' 

#read into dataframe in workspace from bucket
A1c_df <- read_bq_export_from_workspace_bucket(A1c_path)

#check names and ids of measurements #sometimes something unexpected or inconsistent is in the raw data
group_by(A1c_df,measurement_concept_id,standard_concept_name) %>% summarize(n=n())

#if undesired measurement types are present, filter
#for A1c, want to KEEP only concept ids 3004410 and 3005673
A1c_df <- A1c_df[(A1c_df$measurement_concept_id==3004410 | 
                  A1c_df$measurement_concept_id==3005673
                    ), ]

#check units
#this will require manual review, as there are a lot of inconsistent and inconsistently-named units which differ from trait to trait
group_by(A1c_df,unit_source_value) %>% summarize(n=n())

#keep only consistent units
#For A1c: "%", "%{HbA1c}", "%{ofHGB}"
A1c_df <- A1c_df %>% filter(unit_source_value == '%' | unit_source_value == '%{HbA1c}' | unit_source_value == '%{ofHGB}')

#exclude non-physiologic values
#for A1c, this was defined as <1% and >30%
A1c_df <- A1c_df %>% filter(value_as_number>1 & value_as_number<30 )

#keep just the earliest value for each subject
A1c_df$measurement_datetime<-as.Date(A1c_df$measurement_datetime)
A1c_df <- A1c_df %>% group_by(person_id) %>% filter(measurement_datetime == min(measurement_datetime))

#some subjects may still have multiple measurements. Keep just one reading per subject
A1c_df <- distinct(A1c_df,person_id,.keep_all=T)

#view resulting distribution
hist(A1c_df$value_as_number,breaks=100)

#keep just "person_id" and "value_as_number" columns, and rename columns for consistency with UKB dataset
A1c_df <- select(A1c_df,person_id,value_as_number)
colnames(A1c_df)<-c('eid',UKB_code)  #eid in UKB == person_id in AoU


#################################################################################
# CLEAN AND FORMAT CATEGORICAL TRAITS

#Example: ICD codes

#Define path to raw data
ICD_path <- paste(your_bucket,'/data/your_ICD_file_here.csv')


# Read the data directly from Cloud Storage into memory.
# This assumes the data is in .csv format, with person_id column    ####
ICD_df <- read_bq_export_from_workspace_bucket(ICD_path)
ICD_df <- rename(ICD_df,'eid'='person_id') #eid in UKB == person_id in AoU

#check which codes are present and how many of each
group_by(ICD_df, source_concept_code) %>% summarize(n=n())


#Per ICD parent code
#Example: Z85
#define UK Biobank field id that corresponds to your parent code (Z85 in this example)
UKB_code <- 'code_41270_Z850' 

#some entries from AoU may be child codes and will need to be mapped to the single UKB field id 
#Z85 child codes include Z85.040, Z85.05, etc.
#create column for presence or absence of any Z85 code (incl child codes)
ICD_df$code_Z850<-ICD_df$source_concept_code %>% str_detect("Z85")

#make df with just those subjects with Z85 codes, then get rid of duplicate subjects
#add column so subjects with this code marked with '1' 
#when this is joined to other trait/code dataframes, subjects without this code will have NA
Z850_df <- ICD_df[ICD_df$code_Z850=="TRUE",] %>% distinct(eid) %>% mutate({{UKB_code}} := 1)


#################################################################################
# MAKE SINGLE PHENOTYPE DATAFRAME
#The cleaned dataframes for each trait are joined into a single dataframe on "eid" (the person_id)

#Example: join cleaned Z850 and A1c dataframes
temp_Z850_A1c <- full_join(Z850_df, A1c_df,by='eid')

#Join all cleaned trait dataframes (categorical and quantitative) on eid to make final phenotype dataframe

#Final phenotype dataframe will be of the format:
#one row per eid (individual) with a column for each trait
#if an individual doesn't have a value for that trait, field will contain "NA"

#Final pheno dataframe: "all_pheno_df"

#Save to your bucket:
destination_filename <- 'all_pheno_df.csv'
write_excel_csv(all_pheno_df, destination_filename)
system(paste0("gsutil cp ./", destination_filename, " ", your_bucket, "/data/"), intern=T)



