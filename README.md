# Autoflow

![](https://raw.githubusercontent.com/NCBI-Hackathons/Autopipeline/master/images/logo.png)

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.30.0-brightgreen.svg)](https://www.nextflow.io/)

A fast, easy way to present complex bioinformatics pipelines to biologists, built on top of [ATACFlow](https://github.com/NCBI-Hackathons/ATACFlow) and [sequence_handling](https://github.com/MorrellLAB/sequence_handling). This pipeline allows users with basic Linux command line knowledge to quickly customize sequence processing pipelines for their study. This pipeline provides support for using multiple popular tools at each stage of the pipeline.

## Intro
There is currently a plethora of tools that process genetic data at various stages of the workflow. With new tools and scripts popping up everyday, navigating different options for file cleaning, processing, and conversion can be difficult. Amalgamating easy-to-use tools in a centralized, dockerized application can help facilitate reproducibility by alleviating some of the burden on bioinformatists. Our application helps researchers move away from producing countless scripts for data load and processing, automating the entire workflow process by allowing users to choose from a suite of pre-defined functions and tools. Autoflow lets users choose their own adventure after an initial on the fly variant lookup.

### What's the Problem
The problem is that there are too many pipelines to ingest and run data. New pipelines appear every time a lab wants to use a particular set of tools. This prevents pipeline standardization and takes too much time to customize. Often times it is faster for experienced bioinformatists to write new pipelines because previous pipelines may be difficult to undersand and apply.

### How our pipeline addresses the problem
Our pipeline allows users to “build” their own pipelines by selecting the pieces and tools they want to include. Our pipeline goes from Fastq files to variant calls. At each stage in the process our pipeline supports multiple tools so users can pick the most appropriate tool for their organism of study. The pipeline uses NextFlow and Singularity as back-end so it follows a standard syntax and can be modified and re-used easily. This allows biologists with limited computational training to build and customize their own sequencng processing pipelines for their own studies.

### Get in the flow, autoflow!

![](https://raw.githubusercontent.com/NCBI-Hackathons/Autopipeline/master/images/Flowchart_v2.png)

# Installation



# To Do
- [ ] Check if marker exists in dataset before running entire pipeline - **Chaochih**
- [x] Get toy dataset - **Dina**
- [x] Commands to download data from multiple data repositories - **Dina**
    - [x] Wrap commands in script - **Chaochih**
- [ ] Script to parse FastQC output and output warning messages where researcher needs to make decision - **Rob**
- [ ] Figure out how to add additional components to Nextflow pipeline - **Anney**
    - [ ] Including sourcing bash scripts and using functions within bash scripts - **Anney** (**Chaochih** can help direct which functions are necessary to use)
- [x] Install Java to run Nextflow on Google server - **Anney**
- [ ] Adapter trimming script, provide multiple tools as options - **Chaochih**
    - [x] Scythe
    - [x] Trimmomatic
    - [ ] Cutadapt
- [ ] Read mapping, provide multiple tools as options - **Chaochih**
    - [x] BWA
    - [ ] Minimap2
    - [ ] Stampy
- [ ] SAM to BAM/CRAM with Samtools - **TBD**
- [ ] Variant calling - **TBD**
    - [x] GATK
    - [x] Freebayes
    - [ ] SAMTools
    - [ ] Platypus
- [ ] Jupyter notebook - **Devante**
- [ ] GitHub Documentation - **Devante** & **Chaochih**
- [ ] Dockerize pipeline - **TBD**
    - [ ] Do a clean install
- [ ] Manuscript - **Dina**
- [ ] Add support for downstream tools - **TBD**
- [ ] GUI interface to run entire analysis based on parts selected by user (if time) - **Rob**
