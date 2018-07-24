#!/usr/bin/env nextflow

params.config = "$baseDir/config/envVars.sh"
process sourceConfig{
	"""
	source params.config
	"""
}

process QC {
	Main_Quality_Assessment_FastQC 
}