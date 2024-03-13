####################################################################################
#VARIANT ANNOTATION TABLE
####################################################################################
#import and filter the variant annotation table (VAT)
#this VAT is available through AoU and contains annotations including variant consequence (e.g. missense, nonsense), transcript id
#the raw VAT is quite big so it's helpful to filter it for just desired genes, variants before further use
#this also formats and keys variant ids for consistency with WGS/WES matrix tables

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

####################################################################################
#GET GENE LIST
##############
#import table with ensembl ids of your genes of interest
name_of_file_in_bucket = 'gene_df.csv'
os.system(f"gsutil cp '{bucket}/data/{name_of_file_in_bucket}' .")
print(f'[INFO] {name_of_file_in_bucket} is successfully downloaded into your working space')
#read into pandas data frame
gene_df = pd.read_csv(name_of_file_in_bucket)
#convert to Hail table, rename and key fields
gene_ht = hl.Table.from_pandas(gene_df)
gene_ht = gene_ht.rename({'Ensembl_Gene_ID':'gene_id'})
gene_ht = gene_ht.key_by('gene_id')


####################################################################################
#FILTER VAT
###########
vat_path ='X'  #path available in AoU documentation 
vat_table = hl.import_table(vat_path,force=True, delimiter="\t", force_bgz=True, types = {"position":"tint"})

#keep desired subset of columns:
vat_table = vat_table.select('vid', 'gene_symbol','gene_id','transcript','contig','position','ref_allele','alt_allele',
                             'consequence','aa_change','is_canonical_transcript','gnomad_all_af')

#filter for missense variants in canonical transcripts
vat_table = vat_table.filter((vat_table["is_canonical_transcript"]=="true") &
                                          (vat_table.consequence.contains("missense")))
vat_table = vat_table.key_by('gene_id')

#join filtered vat to gene ht on gene ids 
#to get only and all rows in vat that have a gene id in gene_ht
test = gene_ht.join(vat_table, how='left')
#make locus and allele keys from variant ids
test = test.key_by()
test = test.annotate(vid38=test.vid.replace('-',':'))
test = test.annotate(vid38 = hl.str('').join(["chr",test.vid38]))
test = test.key_by(**hl.parse_variant(test.vid38))
#drop unneeded fields
test = test.drop('is_canonical_transcript','vid38','chr','start','end','Gene_Symbol','transcript_id')
#so this is a table of all of the missense variants in the canonical transcripts present in AoU for the genes specified in gene_df

#rename and save to bucket
vat_missense_canonical_99genes = test
vat_missense_canonical_99genes.write(f'{bucket}/data/vat_canonical_missense_99genes.ht')

