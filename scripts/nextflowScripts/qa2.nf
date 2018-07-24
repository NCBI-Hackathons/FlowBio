#!/usr/bin/env nextflow

params.sampleList = "$baseDir/test_data/sample.list"
params.config = "$baseDir/config/envVars.sh"
params.script = "$baseDir/bin/Quality_Assessment.sh"
params.outdir = "$baseDir/results"
params.projectName = "QA"

process sourceFiles{
	"""
	source params.config
	source params.script
	"""
}

process QC {
	Main_Quality_Assessment_FastQC $params.sampleList $params.outdir $params.projectName ${TARGET}
}

workflow.onComplete { 
	println ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir\n" : "Oops .. something went wrong" )
}