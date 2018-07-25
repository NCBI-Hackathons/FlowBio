#!/bin/bash

#   This script uses user defined tools to do read mapping

set -o pipefail

#   What are the dependencies for Read_Mapping?
declare -a Read_Mapping_Dependencies=(bwa parallel)

#   A function to make the output directory
function make_out_dir() {
    local out_dir="$1"
    #   Check if out directory exists, if not make it
    mkdir -p "${out_dir}/Read_Mapping"
}

export -f make_out_dir

#   A function to parse BWA settings
function ParseBWASettings() {
	local POSITIONALS='' # Create a string of positional arguments for BWA mem
	if [[ "${RESCUE}" == true ]]; then POSITIONALS="${POSITIONALS}"'-P '; fi # Add paired to the positionals
	if [[ "${INTERLEAVED}" == true ]]; then POSITIONALS="${POSITIONALS}"'-p '; fi # Add interleaved to the positionals
	if [[ "${SECONDARY}" == true ]]; then POSITIONALS='-a'; else SECONDARY=''; fi # Add secondary to the positionals
	if [[ "${APPEND}" == true ]]; then POSITIONALS="${POSITIONALS}"'-C '; fi # Add append to the positionals
	if [[ "${HARD}" == true ]]; then POSITIONALS="${POSITIONALS}"'-H '; fi # Add hard to the positionals
	if [[ "${SPLIT}" == true ]]; then POSITIONALS="${POSITIONALS}"'-M '; fi # Add split to the positionals
	if [[ "${VERBOSITY}" == 'disabled' ]]; then VERBOSITY=0; elif [[ "${VERBOSITY}" == 'errors' ]]; then VERBOSITY=1; elif [[ "${VERBOSITY}" == 'warnings' ]]; then VERBOSITY=2; elif [[ "${VERBOSITY}" == 'all' ]]; then VERBOSITY=3; elif [[ "${VERBOSITY}" == 'debug' ]]; then VERBOSITY=4; else echo "Failed to recognize verbosity level, exiting..."; exit 1; fi # Set the verbosity level
	MEM_SETTINGS=$(echo "-t ${THREADS} -k ${SEED} -w ${WIDTH} -d ${DROPOFF} -r ${RE_SEED} -A ${MATCH} -B ${MISMATCH} -O ${GAP} -E ${EXTENSION} -L ${CLIP} -U ${UNPAIRED} -T ${RM_THRESHOLD} -v ${VERBOSITY} ${POSITIONALS}") # Assemble our settings
    echo "${MEM_SETTINGS}" # Return our settings
}

#   Export the function
export -f ParseBWASettings

#   A function to see if our referenced FASTA is indexed
function checkIndex() {
    local reference="$1" # What is our reference FASTA file?
    if ! [[ -f "${reference}" ]]; then echo "Cannot find reference genome, exiting..." >&2; exit 1; fi # Make sure it exists
    if ! [[ -r "${reference}" ]]; then echo "Reference genome does not have read permissions, exiting..." >&2; exit 1; fi # Make sure we can read it
    local referenceDirectory=$(dirname "${reference}") # Get the directory for the reference directory
    local referenceName=$(basename "${reference}") # Get the basename of the reference
    if [[ ! $(ls "${referenceDirectory}" | grep "${referenceName}.amb" ) || ! $( ls "${referenceDirectory}" | grep "${referenceName}.ann") || ! $( ls "${referenceDirectory}" | grep "${referenceName}.bwt") || ! $( ls "${referenceDirectory}" | grep "${referenceName}.pac") || ! $( ls "${referenceDirectory}" | grep "${referenceName}.sa") ]]; then return 1; fi # Check that we have all the index files, if we don't then return 1 (not exit 1) so that we can index them
    if [[ ! -r "${referenceDirectory}"/"${referenceName}.amb" || ! -r "${referenceDirectory}"/"${referenceName}.ann" || ! -r "${referenceDirectory}"/"${referenceName}.bwt" || ! -r "${referenceDirectory}"/"${referenceName}.pac" ||! -r "${referenceDirectory}"/"${referenceName}.sa" ]]; then echo "Reference index files do not have read permissions, exiting..." >&2; exit 1; fi # Make sure we can read the index files
}

export -f checkIndex

#   A function to index the FASTA and exit
function indexReference() {
    local reference="$1" # What is our reference FASTA file?
    local referenceDirectory=$(dirname "${reference}") # Get the directory for the reference directory
    if ! [[ -w "${referenceDirectory}" ]]; then echo "Cannot create reference index because you do not have write permissions for ${referenceDirectory}" >&2; exit 1; fi # Make sure we can create the index files
    echo "Indexing reference, will quit upon completion..." >&2
    (set -x; bwa index "${reference}") # Index our reference FASTA file
    echo "Please re-run sequence_handling to map reads" >&2
    exit 10 # Exit the script with a unique exit status
}

export -f indexReference

#   A function to create our read group ID for BWA
function createReadGroupID() {
    local sample="$1" # What is our sample name?
    local project="$2" # What is the name of the project?
    local platform="$3" # What platform did we sequence on?
    local readGroupID="@RG\tID:${sample}\tLB:${project}_${sample}\tPL:${platform}\tSM:${sample}" # Assemble our read group ID
    echo "${readGroupID}" # Return our read group ID
}

export -f createReadGroupID

#   Run read mapping for paired-end samples
function align_bwa_paired() {
    local forward_sample_file="$1" # Where is the forward sample?
    local reverse_sample_file="$2" # Where is the reverse sample?
    local project="$3" # What is the name of our project?
    local platform="$4" # What platform did we sequence on?
    local out_dir="$5"/Read_Mapping # Where is our outdirectory?
    local reference="$6" # Where is our reference FASTA file?
    local mem_settings="$7" # Assemble our settings for BWA mem
    #   Extract sample name from .fq.gz file
    #   Need to modify this part later to take in variations of gzipped fastq file extensions
    #   (i.e. .fastq.gz, .fg.gz, .fq, etc.) or fasta files
    fwd_sample_name=$(basename "${forward_sample_file}" .fq.gz)
    rev_sample_name=$(basename "${reverse_sample_file}" .fq.gz)
    common_name=${fwd_sample_name} # Forward and reverse sample names (may need to tweak this depending on _R1.fastq.gz or _1.fastq.gz)
    #   Assemble our read group ID
    local readGroupID=$(createReadGroupID "${sampleName}" "${project}" "${platform}")
    #   Align sample
    bwa mem "${mem_settings}" -v 2 -R "${readGroupID}" "${reference}" "${fwd_sample_name}" "${rev_sample_name}" > "${out_dir}"/"${common_name}".sam
}

export -f align_bwa_paired

function Main_Read_Mapping_BWA_PE() {
    local forward_list="$1"
    local reverse_list="$2"
    local out_dir="$3"
    #   Parse BWA MEM settings and assemble settings for read alignment
    bwa_mem_settings=$(ParseBWASettings)
    #   Map reads in parallel
    parallel align_bwa_paired {} "${project_name}" "${seq_platform}" "${out_dir}" "${ref}" "${bwa_mem_settings}" :::: "${sample_list}"
}

export -f Main_Read_Mapping_BWA_PE

#   Run read mapping for single-end samples
function align_bwa_singles() {
    local sample_file="$1" # What is the name of our sample?
    local project="$2" # What is the name of our project?
    local platform="$3" # What platform did we sequence on?
    local out_dir="$4"/Read_Mapping # Out dir "Read_Mapping" stores outputs
    local reference="$5" # Where is our reference FASTA file?
    local mem_settings="$6" # Assemble our settings for BWA mem
    #   Extract sample name from .fq.gz file
    #   Need to modify this part later to take in variations of gzipped fastq file extensions
    #   (i.e. .fastq.gz, .fg.gz, .fq, etc.) or fasta files
    sample_name=$(basename "${sample_file}" .fq.gz)
    #   Assemble our read group ID
    local readGroupID=$(createReadGroupID "${sample_name}" "${project}" "${platform}")
    #   Align sample
    bwa mem "${mem_settings}" -v 2 -R "${readGroupID}" "${reference}" "${sampleFile}" > "${out_dir}/${sample_name}.sam"
}

export -f align_bwa_singles

#   Driver function mapping single end reads in parallel with BWA MEM
function Main_Read_Mapping_BWA_SE() {
    local sample_list="$1" # list of samples
    local project_name="$2" # project name used to name summary stats files
    local seq_platform="$3"
    local out_dir="$4"
    local ref="$5"
    #   Parse BWA MEM settings and assemble settings for read alignment
    bwa_mem_settings=$(ParseBWASettings)
    #   Map reads in parallel
    parallel align_bwa_singles {} "${project_name}" "${seq_platform}" "${out_dir}" "${ref}" "${bwa_mem_settings}" :::: "${sample_list}"
}

export -f Main_Read_Mapping_BWA_SE

function align_minimap2_full_genome() {
    local ref="$1"
    local reads="$2"
    local out_dir="$3"
    local preset="$4"
    #   Sample name taken from full name of gzipped FASTA file
    sample_name=$(basename "${reads}" .fa.gz)
    #   Full genome alignment using minimap2
    #   asm10 is one of Minimap2 presets, change depending on organism population diversity
    minimap2 -aLx "${preset}" "${ref}" "${reads}" > "${out_dir}"/"${sample_name}"_asm10.sam
}

export -f align_minimap2_full_genome

function Main_Read_Mapping_Minimap2_FG() {
    local ref="$1"
    local reads_list="$2"
    local out_dir="$3"
    #   Read map our samples using Minimap2 full genome mode
    parallel align_minimap2_full_genome "${ref}" {} "${out_dir}" :::: "${reads_list}"
}

export -f Main_Read_Mapping_Minimap2_FG
