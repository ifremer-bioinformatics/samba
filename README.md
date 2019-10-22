# SAMBA
## Standardized and Automated MetaBarcoding Analysis workflow using Nextflow

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

* q2\_manifest : tabular file with sample name and path to corresponding R1 and R2 fastq.gz files [required]

Mandatory columns are listed below :

| sample-id | forward-absolute-filepath | reverse-absolute-filepath |
| :-------: | :-----------------------: | :-----------------------: |
| sample1   | /path/to/sample1-R1.fastq.gz | /path/to/sample1-R2.fastq.gz |
| sample2   | /path/to/sample2-R1.fastq.gz | /path/to/sample2-R2.fastq.gz |
| sample3   | /path/to/sample3-R1.fastq.gz | /path/to/sample3-R2.fastq.gz |

* q2\_metadata : tabular file describing samples metadata (prefer to use "\_" for long variable names) [required]

Mandatory columns are sample-id and barcode. For your metadata colums, prefer to use "\_" to name your variables :

| sample-id | barcode | your\_metadata\_1 | your\_metadata\_2 |
| :-------: | :-----: | :-------------: | :-------------: | 
| sample1   | ATTAC   | metadata1       | A               |
| sample2   | ACTGA   | metadata1       | B               |
| sample3   | CTTCA   | metadata2       | B               |

* [inasv\_table] : ASV table to use as input if running only statistics steps [optional]

### Workflow parameters

* nextflow.config : general configs for Nextflow (set TMPDIR, Workdir...). Check this file and adapt to your environment !
* config/params.config : workflow workdir definition, processes (tasks) parameters and activation. Check this file and adapt to your data !
* config/resources.config : scheduler resources to attribute to each process. Check this file and adapt to your scheduler !
* config/report.config : nextflow automatic reports parameters 
* SAMBA.nf : each step is described within its command line

### Running the workflow  

* RunSAMBA.sh : Script to run Nextflow command. Check this file and adapt to activate nextflow environment !

### How to get this worflow
```
#connect to datarmor (for Ifremer users)
ssh datarmor
#move to the directory in which you want to get the workflow
cd /TO/YOUR_WORKING_DIRECTORY
## For all users :
#get a local copy of the workflow in the directory SAMBA-nextflow
git clone https://gitlab.ifremer.fr/bioinfo/SAMBA-nextflow
cd SAMBA-nextflow
```
### How to run
Don't forget to modify nextflow config files before running the workflow (see nextflow.config and config directory).

Note : This workflow is design to run on a PBS pro cluster (See resources.config to change cluster options)

```bash
# RUN FROM SCRATCH
./RunSAMBA.sh

# RUN RESUME (when a task has failed or if you run steps separately)
./RunSAMBA.sh -resume
```

### References 

References databases (SILVA 132) are available on our [FTP](ftp://ftp.ifremer.fr/ifremer/dataref/bioinfo/sebimer/sequence-set/qiime2/2019.07/)
Training dataset used from [Qiime2 Tutorial] (https://docs.qiime2.org/2019.7/tutorials/atacama-soils/), [associated publication](https://msystems.asm.org/content/2/3/e00195-16)
