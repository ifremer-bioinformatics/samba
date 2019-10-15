# nextmb
## metabarcoding analysis workflow using nextflow

This workflow uses Qiime 2 (import/cutadapt/dada2/R/phyloseq/vegan) to analyse paired-end metabarcoding data

### Input files 

* q2_manifest : tabular file with sample name and path to corresponding R1 and R2 fastq.gz files
* q2_metadata : tabular file describing samples metadata
* [inasv_table] : ASV table to use as input if running only statistics steps

### Workflow parameters

* config/params.config : workflow workdir definition, processes (tasks) inputs. Check this file and adapt to your data
* config/resources.config : scheduler resources to attribute to each process. Check this file and adapt to your scheduler
* config/report.config : nextflow automatic reports parameters 
* MB.nf : each step is described within its command line


### How to get this worflow
```
#connect to datarmor (for Ifremer users)
ssh datarmor
#move to the directory in which you want to get the workflow
cd /TO/YOUR_WORKING_DIRECTORY
## For all users :
#get a local copy of the workflow in the directory nextmb
git clone https://gitlab.ifremer.fr/bioinfo/nextmb
cd nextmb
```
### How to run
Don't forget to complete nextflow config files before running the workflow (see nextflow.config and config directory)
Note : This workflow is design to run on a grid cluster

```bash
# RUN FROM SCRATCH
./RunNextMB.sh

# RUN RESUME (when a task has failed or if you run steps separately)
./RunNextMB.sh -resume
```
