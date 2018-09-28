#!/usr/bin/env nextflow

Channel
    .fromPath(params.raw_sample_list)
    .splitText(by: 1)
    .set{ sample_paths }

process Quality_Assessment_NanoPlot {
    executor = 'local'
    beforeScript "source ${DEPENDENCIES_FILE}"
    publishDir "${OUTPUT_DIRECTORY}/Quality_Assessment_NanoPlot", mode: 'move', overwrite: false

    input:
    val sample from sample_paths

    output:
    stdout result
    file 'QANP.txt'

    """
    Quality_Assessment_NanoPlot.sh ${sample}
    """
}

result.subscribe { println it }