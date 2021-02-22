#!/bin/bash

# s01_split_multialelic_sites.sh
# Alexey Larionov, 09Feb2021

# Stop at runtime errors
set -e

# Start message
echo "Splitting multiallelic sites"
echo "Started: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# File names
working_folder="/Users/alexey/Documents/mg/s2021/RMSzhang_SJRHB020"
cd "${working_folder}"
source_vcf="${working_folder}/zhang_hg38.bwa.QC.vcf"
split_vcf="${working_folder}/zhang_hg38.bwa.QC.split_MA.vcf"
split_log="${working_folder}/zhang_hg38.bwa.QC.split_MA.log"

# Progress report
echo "Source vcf"
echo "${source_vcf}"
echo ""
echo "Split vcf"
echo "${split_vcf}"
echo ""

# Count variants
echo "Variant counts before splitting"
echo ""
bcftools +counts "${source_vcf}"
echo ""

# Split multiallelic sites
echo "Splitting multiallelic sites ..."

bcftools norm "${source_vcf}" \
--multiallelics -both \
--do-not-normalize \
--output "${split_vcf}" \
&> "${split_log}"

# Count variants
echo ""
echo "Variant counts after splitting"
echo ""
bcftools +counts "${split_vcf}"
echo ""

# Completion mesage
echo "Done: $(date +%d%b%Y_%H:%M:%S)"
echo ""
