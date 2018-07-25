#!/usr/bin/env nextflow

//params.sampleList = "$baseDir/test_data/sample.list"
//params.config = "$baseDir/bin/envVars.sh"
//params.script = "/home/chea/Autopipeline/handlers/Quality_Assessment.sh"
//params.pipeline = "Quality_Assessment"
//params.OUT_DIR = "$baseDir/results"
//params.projectName = "QA"
//params.size = "NA"
//params.sraList = ['SRR5204807','SRR5204808']

// Validate inputs

if (params.PIPELINES){
	//pipeline_list = params.PIPELINES.tokenize(",")

	pipelines = Channel
				.from(params.PIPELINES)
				
}
if ( params.sras ){
	sra_ids_list = params.sras.tokenize(",")
	//if (!params.read_type) exit 1, "Please submit submit your read type of data (paired end - PE or single end - SE)"
} else { 
	Channel.empty().set {sra_ids_list } 
}


process Data_Fetcher{
	publishDir "${params.OUT_DIR}/sra_files/", mode: 'copy'
    tag "reads: ${sra_id}"
	//beforeScript "source Data_Fetcher.sh"
	input:
	val (sra_id) from sra_ids_list

	output:
	set val(sra_id), file("*.fastq") into sra_read_files

	script:
	"""
	source Data_Fetcher.sh
	echo "Downloading ${sra_id} to ${params.OUT_DIR}"
	Main_Fetch_Data ${sra_id} ${params.OUT_DIR}
	echo "Converting SRA to FASTQ in ${params.OUT_DIR} for ${params.read_type}"
	Main_Sra_to_Fastq ${params.read_type} ${params.OUT_DIR}
	#fastq-dump --split-3 ${sra_id}
	# TEST THIS LATER, SHOULD BE FASTER AND DEFAULTS TO --split-3
	#fasterq-dump ${sra_id}
	"""
}
workflow.onComplete { 
	println ( workflow.success ? "\nDone! Open the following report in your browser --> $params.OUT_DIR\n" : "Oops .. something went wrong" )
}

exit 1
process workflow {
	input:
	val(p) from pipelines
	
	script:
	if (p == "Quality_Assessment"){
		if ( params.SAMPLE_LIST ){
    		sampleList = file(params.SAMPLE_LIST)
    		if( !sampleList.exists() ) exit 1, "sample file not found: ${params.SAMPLE_LIST}"
		}
		beforeScript "source ${p}.sh"
		"""
		echo "Run Main_Quality_Assessment_FastQC"
		#echo "source $params.config"
		#echo "source $params.script"
		echo "Main_Quality_Assessment_FastQC $params.SAMPLE_LIST $params.OUT_DIR $params.PROJECT_NAME $params.TARGET"
		"""
	}else if (p == "SayHi"){
		"""
		echo "Hi"
		"""
	}else if (p == "Data_Fetcher"){
		//publishDir "${params.OUT_DIR}/sra_files/", mode: 'copy'
    	tag "reads: ${sra_id}"
		beforeScript "source Data_Fetcher.sh"
    	input:
    	val (sra_id) from sra_ids_list

    	output:
    	set val(sra_id), file("*.fastq") into sra_read_files

    	script:
    	"""
    	#source Data_Fetcher.sh
    	echo "Downloading ${sra_id} to ${params.OUT_DIR}"
    	Main_Fetch_Data ${sra_id} ${params.OUT_DIR}
    	echo "Converting SRA to FASTQ in ${params.OUT_DIR} for ${params.read_type}"
    	Main_Sra_to_Fastq ${params.read_type} ${params.OUT_DIR}
    	#fastq-dump --split-3 ${sra_id}
    	# TEST THIS LATER, SHOULD BE FASTER AND DEFAULTS TO --split-3
    	#fasterq-dump ${sra_id}
   		"""
	}

}



