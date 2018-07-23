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
