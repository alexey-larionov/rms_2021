#!/bin/bash

# s01_copy_vcf_from_aws.sh
# Alexey Larionov, 13Apr2021

# Intended use
# ./s01_copy_vcf_from_aws.sh &> s01_copy_vcf_from_aws.log

# Stop at runtime errors
set -e

# Start message
echo "Copy VCF from AWS"
date
echo ""

# Files and folders
base_folder="/Users/alexey/Documents/mg/s2021/zhang_tests"

scripts_folder="${base_folder}/scripts/s11_vep_annotated_multisample_vcf"
cd "${scripts_folder}"

data_folder="${base_folder}/data/s11_vep_annotated_multisample_vcf"
mkdir -p "${data_folder}"

aws_host="3.10.118.163"
aws_folder="/home/share/data/s03_annotate"
aws_vcf="zhang_hg38.bwa.QC.MA-flag.MA-split.ID.std-Chr.Reheaded.ClinVar.VEP.split-VEP.vcf.gz"

# Progress report
echo "aws_host ${aws_host}"
echo "aws_folder ${aws_folder}"
echo "aws_vcf ${aws_vcf}"
echo "data folder ${data_folder}"
echo ""

# Copy data
echo "Copying ..."
echo ""
rsync -avhe "ssh -i /Users/alexey/.ssh/test01_2020.pem" \
"ec2-user@${aws_host}:${aws_folder}/${aws_vcf}" \
"${data_folder}/"
echo ""

# Completion mesage
echo "Done"
date
echo ""
