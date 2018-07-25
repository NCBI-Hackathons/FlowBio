#!/bin/bash

#   This script generates summary statistics from a VCF file.

set -o pipefail

#   What are the dependencies for Variant_Analysis?
declare -a Variant_Analysis_Other_Dependencies=(python3 libbz2 liblzma freebayes)

#   A function to run each analysis using FreeBayes
function Main_Variant_Analysis_FreeBayes() {
    local fasta_ref="$1" # What is our fasta reference?
    local bam="$2" # What is our bam?
    local out=Variant_Analysis_FreeBayes/$(basename ${2} .vcf)

    #   Make sure the out directory exists
    mkdir -p Variant_Analysis_FreeBayes

    #run the tool
    freebayes --fasta-reference "${fasta_ref}" "${bam}" > "${out}"
}

#   Export the function
export -f Main_Variant_Analysis_FreeBayes

#   A function to run each analysis using SAMtools
function Main_Variant_Analysis_SAMtools() { 

    #   Make sure the out directory exists
}

#   Export the function
export -f Main_Variant_Analysis_SAMtools


#   A function to run each analysis using Platypus
function Main_Variant_Analysis_Platypus() {


    #   Make sure the out directory exists
    mkdir -p "${out}"

}

#   Export the function
export -f Main_Variant_Analysis_Platypus
