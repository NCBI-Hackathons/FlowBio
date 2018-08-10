# To Do List

### General

- [ ] Reorganize Github repo tree structure - assigned to **Chaochih**
- [x] Get toy dataset - assigned to **Dina**
- [x] Draft of manuscript - assigned to **Dina**
- [ ] GUI interface to run entire analysis based on parts selected by user - **Corey Carter**
- [ ] Generate config file from Google Form - assigned to **Corey Carter**
- [ ] Auto-generate config file for user to fill out - assigned to **Corey Carter**
- [ ] Read the docs page for documentation - assigned to **Skylar Wyant** & **Chaochih**
- [ ] Cleanup intermediate files at end of user specified pipeline - assigned to **Corey Carter**

### Nextflow framework

- [x] Figure out how to add additional components to Nextflow pipeline - assigned to **Anney**
    - [x] Including sourcing bash scripts and using functions within bash scripts - **Anney** (**Chaochih** can help direct which functions are necessary to use)
- [ ] Finish coding up Nextflow framework - assigned to **Chaochih**
- [ ] Support for selecting job scheduler in Nextflow pipeline - assigned to **Chaochih** & **TBD**
    - [ ] Support for task arrays on PBS job scheduler - assigned to **Skylar Wyant**
    - [ ] Look into how Nextflow supported batch schedulers work - assigned to **TBD**

### Bash scripts

These scripts will contain functions for modularity. Each tool will have its own function. In some cases, the tool will have multiple functions for different ways to use the tool.
- [x] Commands to download data from SRA - assigned to **Dina**
    - [x] Wrap commands in bash script - assigned to **Chaochih**

- [ ] Script to parse FastQC output and output warning messages where researcher needs to make decision - **TBD**

- [ ] Adapter trimming script, provide multiple tools as options - assigned to **Chaochih**
    - [x] Scythe
    - [x] Trimmomatic
    - [ ] Cutadapt

- [ ] Read mapping, provide multiple tools as options - assigned to **Chaochih**
    - [x] BWA - assigned to **Chaochih**
    - [ ] Minimap2 - assigned to **Chaochih**
    - [ ] Stampy - assigned to **Chaochih**
    - [ ] Bowtie2 - assigned to **TBD**

- [ ] SAM to BAM/CRAM with Samtools - assigned to **Skylar Wyant** & **Chaochih**

- [ ] Variant calling - assigned to **Skylar Wyant** & **Rob**
    - [ ] GATK - assigned to **Skylar Wyant**
    - [x] Freebayes - assigned to **Rob**
    - [x] SAMTools - assigned to **Rob**
    - [x] Platypus - assigned to **Rob**

- [ ] Check if marker exists in dataset before running entire pipeline (SNP calling on the fly) - assigned to **TBD**

- [ ] Add support for downstream tools - assigned to **TBD**

### Containerization

- [ ] Dockerize pipeline - assigned to **Dina** & **Devante**
    - [ ] Do a clean install - assigned to **TBD**
    - [ ] Auto install tools based on ones selected by User - assigned to **TBD**
- [ ] Add support for Singularity - assigned to **TBD**

### Reproducibility

- [ ] Jupyter notebook - assigned to **Devante**
- [ ] GitHub Readme page - assigned to **Chaochih**
