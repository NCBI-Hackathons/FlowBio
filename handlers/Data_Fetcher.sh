#!/bin/bash

set -e
set -o pipefail

#   What are the dependencies for Data_Fetcher.sh?
declare -a Data_Fetcher_Dependencies=(parallel)

#   This script contains functions to take a list of SRA run accessions (SRRxxxxxx)
#   downloads the fastq files and converts them from .sra to .fastq.gz format

#   Function that downloads reads from NCBI SRA data repository
function fetch_data() {
    local acc="$1" # SRRxxxxxx ID as listed in NCBI SRA RunInfo table
    local out_dir="$2" # full path to output directory
    #   Base url listed on NCBI SRA documentation
    local base_url="ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra"
    #   Strep path from accession name
    acc_name=$(basename ${acc})
    #   Build rest of URL
    #   As of 2017-09-17, this format is used:
    #   /sra/sra-instant/reads/ByRun/sra/{SRR|ERR|DRR}/<first 6 characters of accession>/<accession>/<accession>.sra
    query_url="${base_url}/${acc_name:0:3}/${acc_name:0:6}/${acc_name}/${acc_name}.sra"
    #   Go into output directory
    cd "${out_dir}/sra_files"
    wget "${query_url}"
}

export -f fetch_data

function Main_Fetch_Data() {
    local accession_list="$1"
    local out_dir="$2"
    #   Check if out directory exists, if not make it
    mkdir -p "${out_dir}/sra_files"
    #   Download data in parallel
    
    parallel fetch_data {} "${out_dir}" :::: "${accession_list}"
}

export -f Main_Fetch_Data

function Main_Sra_to_Fastq() {
    local lib_layout="$1" # Do we have paired end (PE) or single end (SE) data?
    local out_dir="$2"
    cd "${out_dir}/sra_files"
    sra_array=($(find "${out_dir}/sra_files" -name "*.sra"))
    #   Check if out directory exists, if not make it
    mkdir -p "${out_dir}/raw_fastq"
    #   If samples are paired end, denoted with "PE" as input
    if [ ${lib_layout} == "PE" ]
    then
        #   Split files when converting from .sra to .fastq.gz
        parallel 'fastq-dump --split-files -F --gzip --outdir ${out_dir}/raw_fastq {}' ::: "${sra_array[@]}"
    #   Else, if samples are single end, denoted with "SE" as input
    elif [ ${lib_layout} == "SE" ]
    then
        #   Do not split files when converting from .sra to .fastq.gz
        parallel 'fastq-dump -F --gzip --outdir ${out_dir}/raw_fastq {}' ::: "${sra_array[@]}"
    else
        echo "Please provide valid input for library layout.\n
        Valid input includes: 1) PE and 2) SE"
    fi
}

export -f Main_Sra_to_Fastq
