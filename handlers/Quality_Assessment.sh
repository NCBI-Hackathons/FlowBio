#!/bin/sh

#   This script runs FastQC on a series of samples
#   and writes a summary table of the results

set -e
set -o pipefail

#   What are the dependencies for Quality_Assessment?
declare -a Quality_Assessment_Dependencies=(fastqc parallel)

#   A function to unzip and parse FASTQC files
#   Inspired by Paul Hoffman's RNA version of sequence_handling at https://github.com/LappalainenLab/sequence_handling/
function summarizeQC() {
    local zip_file="$1" # The name of the zip file
    local size="$2" # The estimated size of the covered region
    local out_dir="$3" # The out_dir directory
    local project="$4" # The name of the project
    local sample_name="$(basename ${zip_file} _fastqc.zip)" # The name of the sample
    local zip_dir="$(basename ${zip_file} .zip)" # The name of the directory after unzipping
    unzip -q "${zip_file}" -d "$(dirname ${zip_dir})" # Unzip the zip file
    # Get PASS/WARN/FAIL data from the summary.txt file
    local PerBaseSequenceQuality=$(grep "Per base sequence quality" ${zip_dir}/summary.txt | cut -f 1)
    local PerTileSequenceQuality=$(grep "Per tile sequence quality" ${zip_dir}/summary.txt | cut -f 1)
    local PerSequenceQualityScores=$(grep "Per sequence quality scores" ${zip_dir}/summary.txt | cut -f 1)
    local PerBaseSequenceContent=$(grep "Per base sequence content" ${zip_dir}/summary.txt | cut -f 1)
    local PerSequenceGCContent=$(grep "Per sequence GC content" ${zip_dir}/summary.txt | cut -f 1)
    local PerBaseNContent=$(grep "Per base N content" ${zip_dir}/summary.txt | cut -f 1)
    local SequenceLengthDistribution=$(grep "Sequence Length Distribution" ${zip_dir}/summary.txt | cut -f 1)
    local SequenceDuplicationLevels=$(grep "Sequence Duplication Levels" ${zip_dir}/summary.txt | cut -f 1)
    local OverrepresentedSequences=$(grep "Overrepresented sequences" ${zip_dir}/summary.txt | cut -f 1)
    local AdapterContent=$(grep "Adapter Content" ${zip_dir}/summary.txt | cut -f 1)
    local KmerContent=$(grep "Kmer Content" ${zip_dir}/summary.txt | cut -f 1)
    # Get sequence data from the fastqc_data.txt file
    local ReadCount=$(grep "Total Sequences" ${zip_dir}/fastqc_data.txt | cut -f 2)
    local ReadLength=$(grep "Sequence length" ${zip_dir}/fastqc_data.txt | cut -f 2)
    local GC=$(grep "%GC" ${zip_dir}/fastqc_data.txt | cut -f 2)
    local PercentDeduplicated=$(grep "Total Deduplicated Percentage" ${zip_dir}/fastqc_data.txt | cut -f 2)
    local Encoding=$(grep "Encoding" ${zip_dir}/fastqc_data.txt | cut -f 2)
    # If the size is not set to "NA", calculate read depth estimates
    if [[ "${size}" -ne "NA" ]]
    then
        local LongestRead=$(echo ${ReadLength} | cut -d "-" -f 2) # Sometimes the read length is listed as "1-50", but most of the reads are actually 50
        local ReadDepth=$(( ${ReadCount} * ${LongestRead} / ${size} ))
    else
        local ReadDepth="NA"
    fi
    # Write the sequence data to the summary file
    echo -e "${sample_name}\t${Encoding}\t${ReadLength}\t${ReadCount}\t${ReadDepth}\t${GC}\t${PercentDeduplicated}\t${PerBaseSequenceQuality}\t${PerTileSequenceQuality}\t${PerSequenceQualityScores}\t${PerBaseSequenceContent}\t${PerSequenceGCContent}\t${PerBaseNContent}\t${SequenceLengthDistribution}\t${SequenceDuplicationLevels}\t${OverrepresentedSequences}\t${AdapterContent}\t${KmerContent}" >> "${out_dir}/${project}_quality_summary_unfinished.txt"
    rm -rf "${zip_dir}" # Remove the unzipped directory
    mv "${out_dir}/${sample_name}_fastqc.html" "${out_dir}/HTML_Files/" # Move the HTML file for this sample
    mv "${out_dir}/${sample_name}_fastqc.zip" "${out_dir}/Zip_Files/" # Move the zip file for this sample
}

export -f summarizeQC

#   A function to run quality assessment
function Main_Quality_Assessment_FastQC() {
    local sampleList="$1" # What is our list of samples?
    local out_dir="$2"/Quality_Assessment # Where are we storing our results?
    local project="$3" # What do we call our results?
    local size="$4" # What is the size of the covered region?
    mkdir -p "${out_dir}/HTML_Files" "${out_dir}/Zip_Files" # Make our out_dirput directories
    cat "${sampleList}" | parallel "fastqc --out_dir ${out_dir} {}" # Run FastQC in parallel
    # Make a list of all the zip files
    local zipList=$(find "${out_dir}" -name "*.zip" | sort)
    # Add the header to the quality summary file
    echo -e "Sample name\tEncoding\tRead length\tNumber of reads\tRead depth\t%GC\tDeduplicated percentage\tPer base sequence quality\tPer tile sequence quality\tPer sequence quality scores\tPer base sequence content\tPer sequence GC content\tPer base N content\tSequence length distribution\tSequence duplication levels\tOverrepresented sequences\tAdapter content\tKmer content" > "${out_dir}/${project}_quality_summary_unfinished.txt"
    # Calculate stats and add a row to the summary file for each sample
    parallel -v summarizeQC {} "${size}" "${out_dir}" "${project}" :::: "${zipList}"
    # Add the header to a new file to contain the sorted list
    echo -e "Sample name\tEncoding\tRead length\tNumber of reads\tRead depth\t%GC\tDeduplicated percentage\tPer base sequence quality\tPer tile sequence quality\tPer sequence quality scores\tPer base sequence content\tPer sequence GC content\tPer base N content\tSequence length distribution\tSequence duplication levels\tOverrepresented sequences\tAdapter content\tKmer content" > "${out_dir}/${project}_quality_summary.txt"
    # Sort the summary file based on sample name
    tail -n +2 "${out_dir}/${project}_quality_summary_unfinished.txt" | sort >> "${out_dir}/${project}_quality_summary.txt"
    # Remove the unsorted file
    rm "${out_dir}/${project}_quality_summary_unfinished.txt"
}

export -f Main_Quality_Assessment_FastQC
