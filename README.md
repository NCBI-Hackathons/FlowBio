# Autopipeline
A fast, easy way to present complex bioinformatics pipelines to biologists, built on top of ATACFlow (https://github.com/NCBI-Hackathons/ATACFlow).

# Outline
**Goal**: Build a pipeline packager and allow users to “build” their own pipelines. Biologists will select tools they want to include. Make it easy for biologists to take complex bioinformatics workflows and run them. Auto pipeline will be in a standard framework (i.e. CWL)

# Features 
* Download data from primary data repositories. 
* Multiple tools for each step in the pipeline (provides flexibility and customizability for users).
* Include warning messages throughout for steps where researchers should interpret outputs.
* Ability to run same command and have pipeline pick up where it left off or switch to run each part one at a time or the entire pipeline all at once. 
* Wrapped in Nextflow. 

![auto_pipeline](https://user-images.githubusercontent.com/29574436/43094530-38ad445e-8e81-11e8-8d79-653be0fcd6b7.png)

