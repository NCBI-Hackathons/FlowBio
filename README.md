# Autopipeline
A fast, easy way to present complex bioinformatics pipelines to biologists, built on top of ATACFlow (https://github.com/NCBI-Hackathons/ATACFlow).

# Outline
**Goal**: Build a pipeline packager and allow users to “build” their own pipelines. Biologists will select tools they want to include. Make it easy for biologists to take complex bioinformatics workflows and run them. Auto pipeline will be in a standard framework (i.e. CWL)

# To Do
- [ ] Check if marker exists in dataset before running entire pipeline - **Chaochih**
- [ ] Get toy dataset - **Dina**
- [ ] Script to download data from multiple data repositories - **Dina**
    - [ ] Figure out how to download data from each data repository - **Dina**
- [ ] Script to parse FastQC output and output warning messages where researcher needs to make decision - **Rob**
- [ ] Figure out how to add additional components to Nextflow pipeline - **Anney**
- [ ] Install Java to run Nextflow on Google server - **Anney**
- [ ] Adapter trimming, provide multiple tools as options - **TBD**
    - [ ] Scythe
    - [ ] Trimmomatic
    - [ ] Cutadapt
- [ ] Read mapping, provide multiple tools as options - **TBD**
    - [ ] BWA
    - [ ] Minimap2
    - [ ] Stampy
- [ ] SAM to BAM/CRAM with Samtools - **TBD**
- [ ] Variant calling - **TBD**
    - [ ] GATK
    - [ ] Freebayes
- [ ] Jupyter notebook - **Devante**
- [ ] Dockerize pipeline - **TBD**
- [ ] Add support for downstream tools - **TBD**
- [ ] GUI interface to run entire analysis based on parts selected by user (if time) - **Rob**

![auto_pipeline](https://user-images.githubusercontent.com/29574436/43094530-38ad445e-8e81-11e8-8d79-653be0fcd6b7.png)
