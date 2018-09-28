#!/usr/bin/env bash

set -euo pipefail

# Set variables passed from nextflow
sample="$1"
# # Get the name of the sample without the path in front
# # If it is gzipped, this will strip off the .gz
# base_name=$(basename ${sample} .gz)
# # Get the name of the sample without the file extension
# sample_name=$(basename ${base_name} .fastq)
# # Make the output directory
# mkdir -p "${out}/${sample_name}"
# # Generate some plots
# NanoPlot --outdir "${out}/${sample_name}" --prefix "${sample_name}_" --format pdf --plots dot --plots kde --plots hex --plots pauvre --N50 --title "${project} ${sample_name} Quality Assessment" --fastq_rich "${sample}"

#vcftools
sample="$1"
echo "QANP is running on ${sample}"
cat "${sample}"
touch QANP.txt
