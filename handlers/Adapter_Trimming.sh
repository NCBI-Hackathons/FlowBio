#!/bin/bash

#   This script performs adapter trimming
#   on a series of FASTQ samples using tool of choice

set -e
set -o pipefail

#   What are the dependencies for Adapter_Trimming?
declare -a Adapter_Trimming_Dependencies=(scythe parallel java)

#   A function to make the output directory
function make_out_dir() {
    local out_dir="$1"
    #   Check if out directory exists, if not make it
    mkdir -p "${out_dir}"
}

export -f make_out_dir

#   A function to perform the trimming using Scythe
#   This function takes in single sample, multiple samples are
#   run in parallel with driver function call
#       Adapted from Tom Kono
function trim_adapters_scythe() {
    #   Set the arguments for the trimming
    local sample="$1" # What sample are we working with?
    local out_dir="$2" # out_directory
    local adapters="$3" # Adapter file
    local prior="$4" # Prior
    local platform="$5" # Quality encoding platform
    local forward_naming="$6" # What is the forward naming scheme?
    local reverse_naming="$7" # What is the reverse naming scheme?
    #   Check to see if we have a forward or reverse sample
    if [[ ! -z "${forward_naming}" && $(echo "${sample}" | grep "${forward_naming}") ]] # Is this a forward sample?
    then # If yes
        local name="$(basename ${sample} ${forward_naming})"_Forward # Make the name say forward
    elif [[ ! -z "${reverse_naming}" && $(echo "${sample}"| grep "${reverse_naming}") ]] # Is this the reverse sample?
    then # If yes
        local name="$(basename ${sample} ${reverse_naming})"_Reverse # Make the name say reverse
    else # If this is neither
        local name="$(basename ${sample} | cut -f 1 -d '.')"_Single # Make the name say single
    fi

    #   Is our sample compressed?
    if [[ "$(echo ${sample} | rev | cut -f 1 -d '.' | rev)" == "gz" ]] # Is this gzipped?
    then # If so
        local toTrim="${out_dir}/${name}_PIPE" # Make a name for the pipe which will be passed to the trimmer
        rm -f "${toTrim}" # Remove any pipes with the same name
        mkfifo "${toTrim}" # Make the pipe
        gzip -cd "${sample}" > "${toTrim}" & # Uncompress our sample to the pipe
    elif [[ "$(echo ${sample} | rev | cut -f 1 -d '.' | rev)" == "bz2" ]] # Is this bzipped?
    then # If so
        local toTrim="${out_dir}/${name}_PIPE" # Make a name for the pipe which will be passed to the trimmer
        rm -f "${toTrim}" # Remove any pipes with the same name
        mkfifo "${toTrim}" # Make the pipe
        bzip2 -cd "${sample}" > "${toTrim}" & # Uncompress our sample to the pipe
    else # Otherwise
        local toTrim="${sample}" # Our name will be the sample itself
    fi
    #   Make a named for the trimmed sample
    local trimmed="${out_dir}/${name}_ScytheTrimmed.fastq.gz"
    #   Trim the sample
    scythe -a "${adapters}" -p "${prior}" -q "${platform}" --quiet "${toTrim}" | gzip -c > "${trimmed}"
}

export -f trim_adapters_scythe # export the function

#   A driver function to run the adapter trimming function above in parallel
#   Example Usage: source /autopipeline/handlers/Adapter_Trimming.sh && Main_Adapter_Trimming_Scythe "/autopipeline/data/test_data/list.txt"  "/autopipeline/data"  "TrimmomaticTest" "_R1.fastq.gz" "_R2.fastq.gz" "/autopipeline/data/test_data/accession_list.txt" "0.5" "sanger"
function Main_Adapter_Trimming_Scythe() {
    local raw_samples="$1" # What is our list of samples?
    local out_directory="$2"/Adapter_Trimming # Where are we storing our results?
    local project="$3" # What do we call our results?
    local forward_naming="$4" # What is the extension indicating a forward read?
    local reverse_naming="$5" # What is the extension indicating a reverse read?
    local adapters="$6" # What is our adapter file?
    local prior="$7" # What is Scythe's prior?
    local platform="$8" # What platform did we sequence on?
    #   Make sure the out directory exists
    mkdir -p "${out_directory}"
    if [[ "$?" -ne 0 ]]; then echo "Unbalanced forward and reverse reads" >&2; exit 1; fi # If not an equal amount, exit out with error
    parallel trim_adapters_scythe {} "${out_directory}" "${adapters}" "${prior}" "${platform}" "${forward_naming}" "${reverse_naming}" :::: "${raw_samples}" # Perform the trim
    find "${out_directory}" -type p -exec rm {} \; # Clean up all pipes
    find "${out_directory}" -name "*_ScytheTrimmed.fastq.gz" | sort > "${out_directory}"/"${project}"_trimmed_adapters.txt # Create our list of trimmmed files
}

export -f Main_Adapter_Trimming_Scythe # export the function

#   A function to perform adapter trimming using Trimmomatic
#   Example Usage: source /autopipeline/handlers/Adapter_Trimming.sh && Main_Adapter_Trimming_Trimmomatic "/autopipeline/data/test_data/GRCh38_reference.fa" "/autopipeline/data/test_data/test_final.bam" 
function trim_adapters_trimmomatic() {
    local sample="$1" # What sample are we working with?
    local out_dir="$2" # out_directory
    local adapters="$3" # Adapter file
    local trimmomatic_jar_file="$4" # path to Trimmomatic .jar file
    local lead_qual_below_n="$5" # Remove leading bases w/ quality below n (i.e. 3)
    local trail_qual_below_n="$6" # Remove trailing bases w/ quality below n (i.e. 3)
    local win_size="$7" # scan with specified window size
    local avg_qual_below_n="$8" # cut when avg quality in window size drops below n (i.e. 15)
    local min_len="$9" # drop reads below n bases long (i.e. 36)
    local forward_naming="${10}" # What is the forward naming scheme?
    local reverse_naming="${11}" # What is the reverse naming scheme?
    #   Check to see if we have a forward or reverse sample
    if [[ ! -z "${forward_naming}" && $(echo "${sample}" | grep "${forward_naming}") ]] # Is this a forward sample?
    then # If yes
        local name="$(basename ${sample} ${forward_naming})"_Forward # Make the name say forward
    elif [[ ! -z "${reverse_naming}" && $(echo "${sample}"| grep "${reverse_naming}") ]] # Is this the reverse sample?
    then # If yes
        local name="$(basename ${sample} ${reverse_naming})"_Reverse # Make the name say reverse
    else # If this is neither
        local name="$(basename ${sample} | cut -f 1 -d '.')"_Single # Make the name say single
    fi

    #   Make a named for the trimmed sample
    local trimmed="${out_dir}/${name}_TrimmomaticTrimmed.fastq.gz"

    ########## Still need to add support for paired end samples #############
    ########## and changing trimming parameters #############
    #   Trim the sample
    java -jar "${trimmomatic_jar_file}" SE -phred33 "${sample}" "${out_dir}/${name}_TrimmomaticTrimmed.fastq.gz" ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:"${lead_qual_below_n}" TRAILING:"${trail_qual_below_n}" SLIDINGWINDOW:"${win_size}":"${avg_qual_below_n}" MINLEN:"${min_len}"
}

export -f trim_adapters_trimmomatic

function Main_Adapter_Trimming_Trimmomatic() {
    local raw_samples="$1" # What is our list of samples?
    local out_directory="$2"/Adapter_Trimming # Where are we storing our results?
    local project="$3" # What do we call our results?
    local forward_naming="$4" # What is the extension indicating a forward read?
    local reverse_naming="$5" # What is the extension indicating a reverse read?
    local adapters="$6" # What is our adapter file?
    local prior="$7" # What is Scythe's prior?
    local platform="$8" # What platform did we sequence on?
    if [[ "$?" -ne 0 ]]; then echo "Unbalanced forward and reverse reads" >&2; exit 1; fi # If not an equal amount, exit out with error
    #   Trim samples in parallel
    parallel trim_adapters_trimmomatic {} "${out_directory}" "${adapters}" "${trimmomatic_jar_file}" "${lead_qual_below_n}" "${trail_qual_below_n}" "${win_size}" "${avg_qual_below_n}" "${min_len}" "${forward_naming}" "${reverse_naming}" :::: "${raw_samples}"
    find "${out_directory}" -name "*_TrimmomaticTrimmed.fastq.gz" | sort > "${out_directory}"/"${project}"_trimmed_adapters.txt # Create our list of trimmmed files
}

export -f Main_Adapter_Trimming_Trimmomatic
