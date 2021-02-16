#!/bin/bash

# s05_extract_genes.sh
# Alexey Larionov, 15Feb2021

# Intended use
# ./s01_extract_genes.sh > s01_extract_genes.log

# Stop at runtime errors
set -e

# Start message
echo "Extract gene names"
date
echo ""

# Files and folders
base_folder="/Users/alexey/Documents/mg/s2021/zhang_tests"

working_folder="${base_folder}/s03_explore_gene_names"
cd "${working_folder}"

data_folder="${base_folder}/d03_germline_annotated_files"
samples_file="${data_folder}/zhang.GL.txt"

all_gene_names_file="${working_folder}/all_gene_names_in_zhang_germline_files.txt"
tmp_file="${working_folder}/tmp.txt"

# Progress report
echo "Samples file"
echo "${samples_file}"
echo ""
echo "Gene names"
echo "${all_gene_names_file}"
echo ""

# Read names of germline samples
samples=$(cat "${samples_file}")

# Loop over samples
for sample in $samples
do

  # Get name of the germline annotated file
  source_file="${data_folder}/${sample}_hg38.bwa.QC.scores.filter.txt"

  # Extract gene names
  awk 'BEGIN {FS="\t"} NR>1 {print $5}' "${source_file}" >> "${tmp_file}"

done # Next sample

# Remove duplicates and sort
cat "${tmp_file}" | sort | uniq > "${all_gene_names_file}"

# Count number of genes
echo "Number of detected genes"
cat "${all_gene_names_file}" | wc -l
echo ""

# Remove temporary file
rm "${tmp_file}"

# Completion mesage
echo "Done"
date
echo ""
