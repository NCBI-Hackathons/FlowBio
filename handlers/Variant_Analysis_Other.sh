#!/bin/bash

#   This script generates summary statistics from a VCF file.

set -o pipefail

#   What are the dependencies for Variant_Analysis?
declare -a Variant_Analysis_Other_Dependencies=(python3 libbz2 liblzma freebayes)

#   A function to run each analysis using FreeBayes
#   Example Usage: source /autopipeline/handlers/Variant_Analysis_Other.sh && Main_Variant_Analysis_FreeBayes "/autopipeline/data/test_data/GRCh38_reference.fa" "/autopipeline/data/test_data/test_final.bam" 
function Main_Variant_Analysis_FreeBayes() {
    local fasta_ref="$1" # What is our fasta reference?
    local bam="$2" # What is our bam?
    local out=/autopipeline/data/Variant_Analysis_FreeBayes

    #   Make sure the out directory exists
    mkdir -p "${out}"

    #run the tool
    freebayes --fasta-reference "${fasta_ref}" "${bam}" > "${out}/$(basename ${2} .bam)".vcf
}

#   Export the function
export -f Main_Variant_Analysis_FreeBayes

#   A function to run each analysis using SAMtools
function Main_Variant_Analysis_SAMtools() { 
    local fasta_ref="$1" # What is our fasta reference?
    local bam="$2" # What is our bam?
    local out=/autopipeline/data/Variant_Analysis_SAMtools
    #   Make sure the out directory exists
    mkdir -p "${out}"
    bcftools mpileup -f "${fasta_ref}" "${bam}" | bcftools call -mv -Ob -o "${out}/$(basename ${2} .bam)".bcf
}

#   Export the function
export -f Main_Variant_Analysis_SAMtools


#   A function to run each analysis using Platypus
#   Example Usage: source /autopipeline/handlers/Variant_Analysis_Other.sh && Main_Variant_Analysis_Platypus "/autopipeline/data/test_data/GRCh38_reference.fa" "/autopipeline/data/test_data/test_final.bam"
function Main_Variant_Analysis_Platypus() {
    local fasta_ref="$1" # What is our fasta reference?
    local bam="$2" # What is our bam?
    local out=/autopipeline/data/Variant_Analysis_Platypus

    #   Make sure the out directory exists
    mkdir -p "${out}"

    #run the tool
    python /autopipeline/scripts/Platypus/bin/Platypus.py callVariants --bamFiles="${bam}" --refFile="${fasta_ref}"  --output="${out}/$(basename ${2} .bam)".vcf

}

#   Export the function
export -f Main_Variant_Analysis_Platypus
