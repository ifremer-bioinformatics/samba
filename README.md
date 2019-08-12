# nextmb
##metabarcoding analysis workflow using nextflow

This worflow uses Qiime 2 (import/cutadapt/dada2) to analyse paired-end metabarcoding data

### Input files 
q2_manifest : tabular file with sample name and path to corresponding R1 and R2 fastq.gz files
q2_metadata : tabular file describing samples metadata

* nextflow.config : workflow workdir definition, processes (tasks) inputs and ressources parameters, reports definitions
* MB.nf : each step is describe within its command line

(TO BE COMPLETED)
