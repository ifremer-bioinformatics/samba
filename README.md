# nextmb
## metabarcoding analysis workflow using nextflow

This workflow uses Qiime 2 (import/cutadapt/dada2) to analyse paired-end metabarcoding data

### Input files 

* q2_manifest : tabular file with sample name and path to corresponding R1 and R2 fastq.gz files
* q2_metadata : tabular file describing samples metadata

### Workflow parameters

* nextflow.config : workflow workdir definition, processes (tasks) inputs and resources parameters, reports definitions (Check TO_BE_SET parameters and replace with your own values)
* MB.nf : each step is describe within its command line


### How to get this worflow
```
#connect to datarmor
ssh datarmor
#move to the directory in which you want to get the workflow
cd /TO/YOUR_WORKING_DIRECTORY
#get a local copy of the workflow in the directory nextmb
git clone https://gitlab.ifremer.fr/bioinfo/nextmb
cd nextmb
```
### How to run
Don't forget to complete nextflow config files before running the workflow (see nextflow.config and config directory)

```bash
# RUN FROM SCRATCH
./RunNextMB.sh

# RUN RESUME (when a task has failed or if you run steps separately)
./RunNextMB.sh -resume
```
