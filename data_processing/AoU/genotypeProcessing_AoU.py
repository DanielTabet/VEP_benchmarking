####################################################################################
#This code processes genetic data from All of Us (AoU)
#Uses the AoU whole exome callset (which is a convenience subset of their whole genome callset), 
#This does QC and then generates a csv for each gene of interest with variant and carrier pairs
#Used in Jupyter Notebooks in All of Us workbench (Python kernel)

####################################################################################
#SETUP
######

import os
import subprocess
import pandas as pd

#import Hail, initialize Spark, and set the default reference to GRCh38.
import hail as hl
hl.init(default_reference="GRCh38")

#save the bucket path as a variable
bucket = os.getenv("WORKSPACE_BUCKET")
#save the genomic data bucket location
genomic_location = os.getenv("CDR_STORAGE_PATH")
#save the path to the WGS Hail MatrixTable
mt_wes_path = os.getenv("WGS_EXOME_MULTI_HAIL_PATH")

#import your table with gene names and positions from bucket if not already in workspace
gene_df_path = f'{bucket}/data/gene_df.csv'
os.system(f"gsutil cp {gene_df_path} .")
#read into pandas data frame
gene_df = pd.read_csv('gene_df.csv')

#Get the filtered variant annotation table (VAT)
#VAT has info on variant consequence (e.g. missense), transcript
#VAT from AoU is a huge file so it's easier to filter that first and only once
#See VAT_processing.py script to generate this filtered VAT
vat_table = hl.read_table(f'{bucket}/data/vat_canonical_missense_99genes.ht')

#import hail-formatted gnomad dataset
gnomad='gs://gcp-public-data--gnomad/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.sites.ht'
gnomad_ht = hl.read_table(gnomad)

####################################################################################
#Get genetic data from AoU for each gene and run additional QC
##############################################################

#Get whole exome callset (subset of whole genome callset in AoU)
mt_wes = hl.read_matrix_table(mt_wes_path)

#Filter flagged samples
flagged_samples = f'{genomic_location}/wgs/short_read/snpindel/aux/qc/flagged_samples.tsv'
sample_to_remove = hl.import_table(flagged_samples, key="s")
mt_wes = mt_wes.anti_join_cols(sample_to_remove)

#Generate QC'd, filtered matrix table for each gene in gene_df
n=gene_df.shape[0] #n rows in gene_df

for i in range(0,n): 
    gene = gene_df.Gene_Symbol[i]
    interval = f'chr{gene_df.chr[i]}:{gene_df.start[i]}-{gene_df.end[i]}'
    #make matrix table for this gene
    mt = hl.filter_intervals(mt_wes,[hl.parse_locus_interval(x,) 
                                     for x in [interval]])    
    #split multi-allelic variants
    mt = hl.split_multi_hts(mt)
    #remove variants with no high-quality genotypes
    #criteria: GQ>=20, DP>=10, and AB>=0.2 for heterozygotes
    mt = mt.filter_rows(~mt.filters.contains('NO_HQ_GENOTYPES'))
    #filter by call rate 
    mt= mt.filter_rows(mt.variant_qc.call_rate > 0.90, keep = True)
    #filter by AoU AF
    mt= mt.filter_rows(mt.info.AF[0] < 0.001, keep = True) 
    #keep just snps
    mt = mt.filter_rows(hl.is_snp(mt.alleles[0],mt.alleles[1]))
    #remove entries with GQ score <20 (Phred)
    mt = mt.filter_entries(mt.GQ > 20)
    #join mt and vat so that only variants in the both tables are kept
    #vat table already filtered for missense variants in canonical transcripts
    mt = mt.semi_join_rows(vat_table)
    #add gnomad annotations
    mt = mt.annotate_rows(gnomad=gnomad_ht[mt.locus, mt.alleles])
    #filter on gnomad AF #keep NA's (AF=0) and AF <0.001
    mt = mt.filter_rows((mt.gnomad.freq[0].AF < 0.001) | (~hl.is_defined(mt.gnomad.freq[0].AF)), keep = True)
    mt = mt.drop('gnomad', 'variant_qc')

    print(f'Final {gene} MT dimensions = {mt.count()}')

    mt.write(f'{bucket}/data/missense_only/{gene}_MAF001.mt',overwrite=True)
    


####################################################################################
#GET VARIANT CARRIER IDS and GENOTYPES
######################################

n=gene_df.shape[0] #n rows in gene_df

for i in range(0,n): 
      gene = gene_df.Gene_Symbol[i]
      Ensembl_Gene_ID = gene_df.Ensembl_Gene_ID[i]
      chr_name = f'chr{gene_df.chr[i]}'
      #read in QC'd matrix table
      mt_gene_path = f'{other_bucket}data/missense_only/{gene}_MAF001.mt'  
      mt_gene = hl.read_matrix_table(mt_gene_path)
      #get rid of non-genotype (GT) entries
      mt_gene = mt_gene.select_entries(mt_gene.GT) 
      #convert from matrix table to table
      ht_gene = mt_gene.entries()
      #keep only subjects that are not homozygous for the reference genotype
      ht_gene = ht_gene.filter(ht_gene.GT.is_non_ref()==True)
      #remove pre-existing keys  
      ht_gene = ht_gene.key_by()
      #generate variant string in desired format
      ht_gene = ht_gene.annotate(contig1 = ht_gene.locus.contig)
      ht_gene = ht_gene.annotate(string=hl.variant_str(ht_gene.locus,ht_gene.alleles))
      ht_gene = ht_gene.transmute(variant_name=ht_gene.string.replace(':','_'))
      #keep just variant name and person_id
      ht_gene =ht_gene .select('variant_name','s')
      #now write to dataframe and add gene-specific columns
      ht_gene_df = ht_gene.to_pandas()
      #columns to add: gene_symbol, ensembl_gene_id, 'chr' (will be the same for all rows per gene)  
      ht_gene_df = ht_gene_df.assign(Gene_symbol = gene, Ensembl_id = Ensembl_Gene_ID, chr = chr_name)
      #rename 's' to 'eid' for consistency with UKB datasets 
      ht_gene_df = ht_gene_df.rename(columns={"s":"eid"})   
      #copy as csv to the bucket
      destination_filename = f"{gene}_variants_filtered.csv"
      test_df.to_csv(destination_filename,sep='\t',header=True,index=False,quoting=None)  
      args = ["gsutil", "cp", f"./{destination_filename}", f"{bucket}/data/"]
      output = subprocess.run(args, capture_output=True)

      # print output from gsutil
      output.stderr  

