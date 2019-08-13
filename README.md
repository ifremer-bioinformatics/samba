# nextmb
## metabarcoding analysis workflow using nextflow

This workflow uses Qiime 2 (import/cutadapt/dada2) to analyse paired-end metabarcoding data

### Input files 

* q2_manifest : tabular file with sample name and path to corresponding R1 and R2 fastq.gz files
* q2_metadata : tabular file describing samples metadata

### Workflow parameters

* nextflow.config : workflow workdir definition, processes (tasks) inputs and ressources parameters, reports definitions (Check TO_BE_SET parameters and replace with your own values)
* MB.nf : each step is describe within its command line

### How to run
```bash
ssh datarmor
bash
. /appli/bioinfo/nextflow/19.07.0/env.sh
cd /TO/YOUR_WORKING_DIRECTORY
git clone https://gitlab.ifremer.fr/bioinfo/nextmb
cd nextmb
nextflow run MB.nf
```
to resume (when a task has failed) :

```bash
nextflow run MB.nf -resume
```
