#!/bin/bash

# Set global variables
DIR_PATH=$1
INPUT_VCF_FILENAME_FILE=$2
OUTPUT_PATH=$3
INSTANCE_TYPE=$4

# Loop through each VCF file in the input file to submit jobs
for VCF_FILENAME in $(cat $INPUT_VCF_FILENAME_FILE)
do
    # Get the prefix of the VCF file
    VCF_FILENAME_PREFIX=$(echo $VCF_FILENAME | cut -d'.' -f1)

    # Submit job to run VEP on each VCF file
    JOB_ID=`dx run swiss-army-knife/4.5.1 -iin="$DIR_PATH/$VCF_FILENAME" -iin='/Tools/docker_images/vep.tar.gz' -icmd="chmod a+rwx /home/dnanexus/out/out && docker load --input vep.tar.gz && docker run -w /tmp -v /home/dnanexus/out/out:/tmp -v /mnt/project:/mnt/project ensemblorg/ensembl-vep:release_105.0 bash -c 'vep -i $VCF_FILENAME_PREFIX.vcf.gz -o $VCF_FILENAME_PREFIX.txt --stats_file $VCF_FILENAME_PREFIX.html  --cache --dir_cache /mnt/project/Caches/vep_data/ --dir_plugins /mnt/project/Caches/vep_data/Plugins --offline --canonical --force_overwrite --af_gnomad --buffer_size 100 --verbose --plugin dbNSFP,/mnt/project/Caches/vep_data/Plugins/dbNSFP4.2a_grch38.gz,Uniprot_acc,hg19_chr,hg19_pos\(1-based\),hg18_chr,hg18_pos\(1-based\),genename,aaref,aaalt,aapos,HGVSp_VEP,HGVSp_ANNOVAR,gnomAD_exomes_AF,gnomAD_genomes_AF,clinvar_hgvs,clinvar_review,clinvar_id,clinvar_clnsig,Ensembl_proteinid,Ensembl_geneid,Ensembl_transcriptid,SIFT_score,Polyphen2_HDIV_score,Polyphen2_HVAR_score,LRT_score,MutationTaster_score,FATHMM_score,PROVEAN_score,VEST4_score,MetaSVM_score,MetaLR_score,M-CAP_score,REVEL_score,MutPred_score,MVP_score,MPC_score,PrimateAI_score,DEOGEN2_score,BayesDel_addAF_score,BayesDel_noAF_score,ClinPred_score,LIST-S2_score,CADD_raw,DANN_score,fathmm-MKL_coding_score,fathmm-XF_coding_score,Eigen-raw_coding,Eigen-PC-raw_coding,GenoCanyon_score,integrated_fitCons_score,GM12878_fitCons_score,H1-hESC_fitCons_score,HUVEC_fitCons_score,GERP++_RS,phyloP100way_vertebrate,phyloP17way_primate,phyloP30way_mammalian,phastCons100way_vertebrate,phastCons17way_primate,phastCons30way_mammalian,SiPhy_29way_pi,SiPhy_29way_logOdds'" --destination="$OUTPUT_PATH" --instance-type=$INSTANCE_TYPE --name="VCF to VEP" -y --brief --priority=low`

    # Print job ID to log file
    echo "Job submitted. VCF file: $VCF_FILENAME; ID: $JOB_ID"

done

echo "Submitted all jobs."