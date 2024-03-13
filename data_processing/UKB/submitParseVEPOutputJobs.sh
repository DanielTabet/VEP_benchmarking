#!/bin/bash

# Set global variables
VCF_PATH=$1
INPUT_VCF_FILENAME_FILE=$2
VEP_PATH=$3
OUTPUT_PATH=$4
INSTANCE_TYPE=$5

# Loop through each VCF file in the input file to submit jobs
for VCF_FILENAME in $(cat $INPUT_VCF_FILENAME_FILE)
do
    # Get the prefix of the VCF file
    VCF_FILENAME_PREFIX=$(echo $VCF_FILENAME | cut -d'.' -f1)

    # Submit job to run VEP on each VCF file
    JOB_ID=`dx run swiss-army-knife/4.5.1 -iin="/Analysis Workflow/parse_vep_output.py" -iin="$VEP_PATH/$VCF_FILENAME_PREFIX.txt" -iin="$VCF_PATH/$VCF_FILENAME" -icmd="pip3 install pandas && python3 parse_vep_output.py -vep $VCF_FILENAME_PREFIX.txt -vcf $VCF_FILENAME" --destination="$OUTPUT_PATH" --instance-type=$INSTANCE_TYPE --name="Parse VEP results" -y --brief --priority=low`

    # Print job ID to log file
    echo "Job submitted. VCF file: $VCF_FILENAME; ID: $JOB_ID"

done

echo "Submitted all jobs."