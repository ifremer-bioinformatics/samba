# nf-core/samba: Usage

## Table of contents

* [Table of contents](#table-of-contents)
* [Introduction](#introduction)
* [Running the pipeline](#running-the-pipeline)
  * [Updating the pipeline](#updating-the-pipeline)
  * [Reproducibility](#reproducibility)
* [Main arguments](#main-arguments)

  Mandatory arguments:
  * [`--input_metadata`](#--input_metadata)
  * [`--input_manifest`](#--input_manifest)
  * [`-profile`](#-profile)

  Generic arguments:
  * [`--singleEnd`](#--singleEnd)

  Other options:
  * [`--outdir`](#--outdir)
  * [`-w/--work-dir`](#-w/--work-dir)
  * [`--email`](#--email)
  * [`--email_on_fail`](#--email_on_fail)
  * [`-name`](#-name)
  * [`--projectName`](#--projectName)

  Data integrity:
  * [`--data_integrity_enable`](#--data_integrity_enable)
  * [`--barcode_filter`](#--barcode_filter)
  * [`--primer_filter`](#--primer_filter)
  
  Raw reads cleaning:
  * [`--primerF`](#--primerF)
  * [`--primerR`](#--primerR)
  * [`--errorRate`](#--errorRate)
  * [`--overlap`](#--overlap)
  
  ASVs inference:
  * [`--FtrimLeft`](#--FtrimLeft)
  * [`--RtrimLeft`](#--RtrimLeft)
  * [`--FtruncLen`](#--FtruncLen)
  * [`--RtruncLen`](#--RtruncLen)
  * [`--FmaxEE`](#--FmaxEE)
  * [`--RmaxEE`](#--RmaxEE)
  * [`--minQ`](#--minQ)
  * [`--chimeras`](#--chimeras)

  Merge ASVs tables:
  * [`--dada2merge`](#--dada2merge)
  * [`--merge_tabledir`](#--merge_tabledir)
  * [`--merge_repseqsdir`](#--merge_repseqsdir)

  Distribution based-clustering: 
  * [`--dbotu3_enable`](#--dbotu3_enable)
  * [`--gen_crit`](#--gen_crit)
  * [`--abund_crit`](#--abund_crit)
  * [`--pval_crit`](#--pval_crit)

  Taxonomic assignation:
  * [`--extract_db`](#--extract_db)
  * [`--seqs_db`](#--seqs_db)
  * [`--taxo_db`](#--taxo_db)
  * [`--database`](#--database)
  * [`--confidence`](#--confidence)

  Decontamination:
  * [`--microDecon_enable`](#--microDecon_enable)
  * [`--control_list`](#--control_list)
  * [`--nb_controls`](#--nb_controls)
  * [`--nb_samples`](#--nb_samples)

  Predict functionnal abundance:
  * [`--picrust2_enable`](#--picrust2_enable)
  * [`--method`](#--method)
  * [`--nsti`](#--nsti)

  Differential abundance testing:
  * [`--ancor_var`](#--ancor_var)

  Statistics:
  * [`--stats_alpha_enable`](#--stats_alpha_enable)
  * [`--stats_beta_enable`](#--stats_beta_enable)
  * [`--stats_desc_comp_enable`](#--stats_desc_comp_enable)

  * [`--kingdom`](#--kingdom)
  * [`--taxa_nb`](#--taxa_nb)
  * [`--alpha_div_group`](#--alpha_div_group)
  * [`--beta_div_var`](#--beta_div_var)
  * [`--sets_analysis_crit`](#--sets_analysis_crit)
  * [`--hc_method`](#--hc_method)
  
  * [`--stats_only`](#--stats_only)
  * [`--inasv_table`](#--inasv_table)
  * [`--innewick`](#--innewick)

  Final analysis report:
  * [`--report_enable`](#--report_enable)

## Introduction

Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline

The most simple command for running the pipeline is as follows:

```bash
nextflow run nf-core/samba -profile conda,test
```

This will launch the pipeline with the `conda` and `test` configuration profiles. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/samba
```

### Reproducibility

It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/samba releases page](https://github.com/nf-core/samba/releases) and find the latest version number - numeric only (eg. `2.0.0`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 2.0.0`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For running a (not necessarily stable) development version of the pipeline you can use the `dev` branch with `-r dev`.

## Main arguments

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments. Note that multiple profiles can be loaded, for example: `-profile docker` - the order of arguments is important!

If `-profile` is not specified at all the pipeline will be run locally and expects all software to be installed and available on the `PATH`.

* `conda`
  * A generic configuration profile to be used with [conda](https://conda.io/docs/)
  * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `docker`
  * A generic configuration profile to be used with [Docker](http://docker.com/)
  * Pulls software from dockerhub: [`nfcore/rnaseq`](http://hub.docker.com/r/nfcore/rnaseq/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
  * Pulls software from DockerHub: [`nfcore/rnaseq`](http://hub.docker.com/r/nfcore/rnaseq/)

Profiles are also available to configure the samba workflow and can be combined with execution profiles listed above.

* `test`
  * A profile with a complete configuration for automated testing
  * Includes training dataset so needs no other parameters
* `custom`
  * A profile to complete according to your data and experiment.

### `--input_metadata`

Path to input file with project samples metadata (csv format). Headers of metadata file must follow the Qiime2 requirements [Qiime2 metadata](https://docs.qiime2.org/2019.10/tutorials/metadata/)

### `--input_manifest`

Path to input file with samples reads files paths (csv format). Headers of manifest file must follow the Qiime2 requirements [Qiime2 manifest](https://docs.qiime2.org/2019.10/tutorials/importing/#manifest-file)

## Generic arguments

### `--singleEnd`

Set to true to specify that the inputs are single-end reads. Default is paired-End.

### `--projectName`

Name of the project being analyzed.

## Data integrity process arguments

This process is optional and checks if input datasets are correctly demultiplexed, if primers ratio is high enough, if metadata file is well-formed and creates a CSV report.

### `--data_integrity_enable`

Data integrity checking step. Set to false to deactivate this step. (default = true)

### `--barcode_filter`

Percentage of sample barcode supposed to be found in raw reads (default : 90).

### `--primer_filter`

Percentage of primers supposed to be found in raw reads (default : 70).

## Raw reads cleaning process parameters

### `--primerF` 

Forward primer (to be used in Cutadapt cleaning step).

### `--primerR`

Reverse primer (to be used in Cutadapt cleaning step).

### `--errorRate`

Cutadapt error rate allowed to match primers (default : 0.1).

### `--overlap`

Cutadapt overlaping length between primer and read (default : 18). 

## ASVs inference process arguments

This process is based on [Qiime2/Dada2](https://docs.qiime2.org/2019.10/plugins/available/dada2/).

### `--FtrimLeft`

The number of nucleotides to remove from the start of each forward read (default : 0 = no trimming).

### `--RtrimLeft`

The number of nucleotides to remove from the start of each reverse read (default : 0 = no trimming).

### `--FtruncLen`

Truncate forward reads after FtruncLen bases. Reads shorter than this are discarded (default : 0 = no trimming).

### `--RtruncLen`

Truncate reverse reads after RtruncLen bases. Reads shorter than this are discarded (default : 0 = no trimming).

### `--FmaxEE`

Forward reads with higher than maxEE "expected errors" will be discarded (default = 2).

### `--RmaxEE`

Reverse with higher than maxEE "expected errors" will be discarded (default = 2).

### `--minQ`

After truncation, reads contain a quality score less than minQ will be discarded (default = 10).

### `--chimeras`

Chimera detection method : default = "consensus". Set to "pooled" if the samples in the sequence table are all pooled together for bimera identification.

## Merge ASVs tables process 

This process is optionnal and based on [Qiime2/feature-table function](https://docs.qiime2.org/2019.10/plugins/available/feature-table/). The workflow can begin at this step if you already have Dada2 ASVs tables that you want to merge to perform the analysis.

### `--dada2merge`

Set to true to merge Dada2 ASVs tables.

### `--merge_tabledir`

Path to the directory containing the ASVs tables to merge (this directory must contain only the ASVs tables to merge).

### `--merge_repseqsdir`

Path to the directory containing the representative sequences to merge (this directory must constain only the representative sequences to merge).

## Distribution based-clustering process arguments

This step, based on [dbotu3](https://github.com/swo/dbotu3) is optional if you don\`t want to cluster your ASVs.

### `--dbotu3_enable`

Distribution based-clustering step. Set to false to deactivate this step. (default = true)

### `--gen_crit`

dbotu3 Genetic criterion (default = 0.1).

### `--abund_crit`

dbotu3 Abundance criterion (default = 10).

### `--pval_crit`

dbotu3 P-value criterion (default = 0.0005).

## Taxonomic assignation process arguments

This process is based on [Qiime2/feature-classifier function](https://docs.qiime2.org/2019.10/tutorials/feature-classifier/).

### `--extract_db`

Set to true to extract specific region from reference database (default = false).

### `--seqs_db`

Path to reference database (required if extract_db = true).

### `--taxo_db`

Path to taxonomic reference database (required if extract_db = true).

### `--database`

Path to preformatted QIIME2 format database (required if extract_db = false).

### `--confidence`

RDP confidence threshold (default = 90).

## Decontamination process arguments

This step is optional and based on [microDecon](https://github.com/donaldtmcknight/microDecon) package.

### `--microDecon_enable`

Sample decontamination step. Set to true to activate this step. (default = false)

### `--control_list`

Comma separated list of control samples (e.g : "sample1,sample4,sample7") (required if microDecon_enable = true).

### `--nb_controls`

Number of controled samples listed (required if microDecon_enable = true).

### `--nb_samples`

Number of samples that are not control samples (required if microDecon_enable = true).

## Predict functionnal abundance process arguments

This step is optional and based on [Qiime2/PICRUSt2](https://github.com/gavinmdouglas/q2-picrust2)

### `--picrust2_enable`

Set to true to enable functionnal prediction step. (default = false)

### `--method`

HSP method of your choice. (default = 'mp' ) The most accurate prediction methode. Faster method: 'pic'.

### `--nsti`

Max nsti value accepted. (default = 2) NSTI cut-off of 2 should eliminate junk sequences.

## Differential abundance testing process arguments

Step based on [Qiime2/Composition ancom](https://docs.qiime2.org/2020.2/plugins/available/composition/ancom/)

### `--ancor_var`

According to your metadata file, select the column name corresponding to the variable to group samples for ANCOM analysis. 

## Statistics processes arguments

Alpha diversity, Beta diversity and Descriptive comparisons statistics can be enabled or disabled.
Statistics steps can also being run alone (without the above bioinformatics steps). See below.

### `--stats_alpha_enable`

Set to false to deactivate Alpha diversity statistics step. (default = true)

### `--stats_beta_enable`

Set to false to deactivate Beta diversity statistics steps. (default = true)

### `--stats_desc_comp_enable`

Set to false to deactivate Descriptive comparisons steps. (default = true)

### `--kingdom`

Kingdom to be displayed in barplots. (default = "Bacteria")

### `--taxa_nb`

Number of taxa to be displayed in barplots.

### `--alpha_div_group`

According to your metadata file, select the column name corresponding to the variable to group samples for Alpha diversity.

### `--beta_div_var`

According to your metadata file, select the column name corresponding to the variable of interest for Beta diversity.

### `--sets_analysis_crit`

According to your metadata file, select the column name corresponding to the variable of interest for Descriptive comparisons graphs.

### `--hc_method`

Hierarchical clustering method (default = 'ward.D2').

## Run the statistics steps only processes arguments

Statistics steps can be run without running previous bioinformatics steps. Above statistics arguments must be set to perform statistics only steps.

### `--stats_only`

Perform only statistical analysis (ASV table and newick tree required). Set to true to activate. (default = false)

### `--inasv_table`

if stats_only is activated, set the path to your own ASV table in tsv format.

### `--innewick`

if stats_only is activated, set the path to your own phylogenetic tree in newick format.

## Final report process arguments

This step is optional and create a HTML report of samba analyses.

### `--report_enable`

Set to false to deactivate report creation. (default = true)


## Job resources

### Automatic resubmission

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests

Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files hosted at [`nf-core/configs`](https://github.com/nf-core/configs/tree/master/conf) for examples.

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition below). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack).

## Other command line parameters

### `--outdir`

The output directory where the results will be saved.

### `-w/--work-dir`

The temporary directory where intermediate data will be saved.

### `--email`

Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits.

### `--email_on_fail`

Same as --email, except only send mail if the workflow is not successful.

### `-name`

Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`

Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override pipeline defaults.

### `--custom_config_version`

Provide git commit id for custom Institutional configs hosted at `nf-core/configs`. This was implemented for reproducibility purposes. Default: `master`.

```bash
## Download and use config file with following git commid id
--custom_config_version d52db660777c4bf36546ddb188ec530c3ada1b96
```

### `--custom_config_base`

If you're running offline, nextflow will not be able to fetch the institutional config files
from the internet. If you don't need them, then this is not a problem. If you do need them,
you should download the files from the repo and tell nextflow where to find them with the
`custom_config_base` option. For example:

```bash
## Download and unzip the config files
cd /path/to/my/configs
wget https://github.com/nf-core/configs/archive/master.zip
unzip master.zip

## Run the pipeline
cd /path/to/my/data
nextflow run /path/to/pipeline/ --custom_config_base /path/to/my/configs/configs-master/
```

> Note that the nf-core/tools helper package has a `download` command to download all required pipeline
> files + singularity containers + institutional configs in one go for you, to make this process easier.

### `--max_memory`

Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'`

### `--max_time`

Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--max_cpus`

Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

### `--plaintext_email`

Set to receive plain-text e-mails instead of HTML formatted.

### `--monochrome_logs`

Set to disable colourful command line output and live life in monochrome.

