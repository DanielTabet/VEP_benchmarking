#############################################################################################
#Generate phenotype file for each gene
#Uses cleaned and formatted phenotype data from 'TraitProcessingAoU'
#These will be input for benchmarking scripts

library(tidyverse)

#get path to Cloud Storage bucket associated with your workspace
your_bucket <- Sys.getenv('WORKSPACE_BUCKET')

geno_pheno <- read.csv('variant-effect-predictor-assessment/common/genes.csv',header=T)

#if all_pheno_df is not already in the workspace
all_pheno_df_path <- paste(your_bucket,'/data/all_pheno_df.csv')
system(paste0("gsutil cp ", all_pheno_df_path, " ."), intern=T)
all_pheno_df <- read.csv('all_pheno_df.csv',header=T)

#for each gene/line in geno/pheno
len=dim(geno_pheno)[1]

for (i in 1:len){
    gene <- geno_pheno[i,1]
    filename <- paste(gene,'_phenotypes_filtered.csv', sep='')

    #get and format list of pheno codes for this gene
    phenos <- geno_pheno[i,3]
    #number of phenotypes for this gene
    n<-str_count(phenos,'x') #x as delimiter
    #now we want each pheno code. 
    test2<-str_replace_all(phenos,'x','code_')
    test3<-str_split_fixed(test2,'\n',n)
    print(paste(gene,test3))
    #Then, subset the columns of the measurement dataframe to just those matching these codes
    test4<-all_pheno_df %>% select(eid,any_of(test3[1:n]))
    #keep only rows with at least one non-na entry, not counting person_id column
    test9 <-test4 %>% filter_at(vars(starts_with('code')),any_vars(!is.na(.)))

    #reformat column names
    temp_names <- names(test9)
    new_names <- sub('code_','',temp_names)
    names(test9) <- new_names
    head(test9,20)
    dim(test9)
    
    #name this gene's dataframe
    assign(paste(gene,'_pheno_df',sep=''),test9)
    print(paste(gene,'_pheno_df has dimensions ',dim(test9)[1],'x',dim(test9)[2], sep=''))
    
    #write to csv and save to bucket
    write_excel_csv(test9, filename)
    system(paste0("gsutil cp ./", filename, " ", my_bucket, "/data/"), intern=T)

    #optional 
    #write to local 'input' folder for use with subsequent benchmarking scripts 
    input_path <- paste("variant-effect-predictor-assessment/input",gene,sep='/')
    #make input directory for this gene
    system(paste('mkdir',input_path), intern=T)
    #write gene's phenotype file to its input directory
    system(paste0("gsutil cp ", filename," ", input_path), intern=T)
}