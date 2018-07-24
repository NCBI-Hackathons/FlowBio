# Autopipeline
A fast, easy way to present complex bioinformatics pipelines to biologists, built on top of ATACFlow (https://github.com/NCBI-Hackathons/ATACFlow).

# Outline
**Goal**: Build a pipeline packager and allow users to “build” their own pipelines. Biologists will select tools they want to include. Make it easy for biologists to take complex bioinformatics workflows and run them. Auto pipeline will be in a standard framework (i.e. CWL)

# To Do
- [ ] Check if marker exists in dataset before running entire pipeline - **Chaochih**
- [x] Get toy dataset - **Dina**
- [ ] Commands to download data from multiple data repositories - **Dina**
    - [ ] Wrap commands in script - **Chaochih**
- [ ] Script to parse FastQC output and output warning messages where researcher needs to make decision - **Rob**
- [ ] Figure out how to add additional components to Nextflow pipeline - **Anney**
    - [ ] Including sourcing bash scripts and using functions within bash scripts - **Anney** (**Chaochih** can help direct which functions are necessary to use)
- [x] Install Java to run Nextflow on Google server - **Anney**
- [ ] Adapter trimming script, provide multiple tools as options - **Chaochih**
    - [x] Scythe
    - [x] Trimmomatic
    - [ ] Cutadapt
- [ ] Read mapping, provide multiple tools as options - **Chaochih**
    - [ ] BWA
    - [ ] Minimap2
    - [ ] Stampy
- [ ] SAM to BAM/CRAM with Samtools - **TBD**
- [ ] Variant calling - **TBD**
    - [ ] GATK
    - [ ] Freebayes
- [ ] Jupyter notebook - **Devante**
- [ ] GitHub Documentation - **Devante** & **Chaochih**
- [ ] Dockerize pipeline - **TBD**
- [ ] Manuscript - **Dina**
- [ ] Add support for downstream tools - **TBD**
- [ ] GUI interface to run entire analysis based on parts selected by user (if time) - **Rob**



![auto_pipeline](https://user-images.githubusercontent.com/29574436/43094530-38ad445e-8e81-11e8-8d79-653be0fcd6b7.png)
