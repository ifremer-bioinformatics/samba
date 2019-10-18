# nextmb
## metabarcoding analysis workflow using nextflow

This workflow will process paired-end metabarcoding data. 

* Data integrity process 
    * Checking pourcentage of correct barcode in R1 and R2 files 
    * Checking number of forward and reverse reads by sample
    * Checking number of forward and reverse primers by sample
    * Checking if the same sequencer identifier is present over the data
* Data import process
    * Import data and create objects for Qiime2 analysis
    * Evaluate data quality
* Trimming process
    * Use Qiime2 cutadapt to remove primers from data
* Sample inference process using Dada2
    * Extract ASVs (Amplicon Sequence Variants) from samples using Dada2
    * Create ASVs counting table
* Taxonomy assignment process
    * Use Qiime2 to assign each ASV according to a reference database using RDP
* Statistical analysis process
    * A step to prepare data for stats will create R Phyloseq object for downstream analysis
    * Alpha diversity boxplots and taxonomic graphs are created using R Phyloseq
    * Beta diversity NMDS plots are generated according to several normalisation processes (no-normalisation, standard rarefaction, deseq2 normalisation, CSS normalisation)
* Reporting
    * A report folder will provide the metabarcoding workflow results
    * A workflow execution synthesis will be generated using Nextflow native DAG, timeline, trace and html report

### Input files 

* q2_manifest : tabular file with sample name and path to corresponding R1 and R2 fastq.gz files [required]
* q2_metadata : tabular file describing samples metadata (prefer to use "\_" for long variable names) [required]
* [inasv_table] : ASV table to use as input if running only statistics steps [optional]

### Workflow parameters

* nextflow.config : general configs for Nextflow (set TMPDIR, Workdir...). Check this file and adapt to your environment !
* config/params.config : workflow workdir definition, processes (tasks) parameters and activation. Check this file and adapt to your data !
* config/resources.config : scheduler resources to attribute to each process. Check this file and adapt to your scheduler !
* config/report.config : nextflow automatic reports parameters 
* MB.nf : each step is described within its command line

### Running script 

* RunNextMB.sh : Script to run Nextflow command. Check this file and adapt to activate nextflow environment !

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
Don't forget to modify nextflow config files before running the workflow (see nextflow.config and config directory).

Note : This workflow is design to run on a PBS pro cluster (See resources.config to change cluster options)

```bash
# RUN FROM SCRATCH
./RunNextMB.sh

# RUN RESUME (when a task has failed or if you run steps separately)
./RunNextMB.sh -resume
```
