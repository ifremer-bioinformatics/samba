# samba: Quickstart guide

## Input file requirements

Input data for SAMBA can be provided in two different formats: either a **single XLS file** (EXCEL 97-2004) with two sheets named "manifest" and "metadata" or **2 tsv files**.

### 1 - The XLS sample file

The XLS file must contains two sheets named "manifest" and "metadata". The two sheets must be formatted as describe below for the two TSV files.
 
### 2 - Manifest file

The manifest file give the full path to input file(s) (FASTQ) for each sample (tsv format -> field separator = tabulation). It's consists of 2 (single-end) or 3 columns (paired-end).  The first line (header) of manifest file must follow the Qiime2 requirements.

Example of manifest file for single-end data:
```
sample-id   absolute-filepath
control     /path/to/data/control_R1_001.fastq.gz
BAQ4166.1   /path/to/data/BAQ4166.1_24_L001_R1_001.fastq.gz
BAQ4166.2   /path/to/data/BAQ4166.2_1_L001_R1_001.fastq.gz
BAQ4166.3   /path/to/data/BAQ4166.3_7_L001_R1_001.fastq.gz
BAQ4697.1   /path/to/data/BAQ4697.1_12_L001_R1_001.fastq.gz
BAQ4697.2   /path/to/data/BAQ4697.2_37_L001_R1_001.fastq.gz
```

And for paired-end data:
```
sample-id   forward-absolute-filepath                       reverse-absolute-filepath
control     /path/to/data/control_L001__R1_001.fastq.gz     /path/to/data/control_L001__R2_001.fastq.gz
BAQ4166.1   /path/to/data/BAQ4166.1_24_L001_R1_001.fastq.gz /path/to/data/BAQ4166.1_24_L001_R2_001.fastq.gz
BAQ4166.2   /path/to/data/BAQ4166.2_1_L001_R1_001.fastq.gz  /path/to/data/BAQ4166.2_1_L001_R2_001.fastq.gz
BAQ4166.3   /path/to/data/BAQ4166.3_7_L001_R1_001.fastq.gz  /path/to/data/BAQ4166.3_7_L001_R2_001.fastq.gz
BAQ4697.1   /path/to/data/BAQ4697.1_12_L001_R1_001.fastq.gz /path/to/data/BAQ4697.1_12_L001_R2_001.fastq.gz
```

As we can see, the input data must be compressed in **gzip format**. Moreover the FASTQ file need to be, at least for data integrity step, encoded in Illumina **Casava 1.8**. In this format, the '@' line (sequence identifier) looks like:
```
@M00176:65:000000000-A41FR:1:2114:9875:23134 1:N:0:CAACTAGA
```
Where each field is describe as follow:
```
Intrument name    Run id  Flowcell id   Lane   Tile  Coordinate (x:y)
             \        |        |         |    /       | 
              _______ __ _______________ _ ____ __________ ______________
              @M00176:65:000000000-A41FR:1:2114:9875:23134 1:N:0:CAACTAGA
                                                          /   \         \
                                                      Pair   Filter   Index seq.
```

### 3 - Metadata file
The metadata file give information about your sequencing parameters (barcode and primers) and biological/environment/etc conditions. The number of descriptive fields is unlimited. As for manifest file, metadata is in tsv format and the header must follow the Qiime2 requirements.

Here, an example for paired-end data:

```
sampleid    barcode         PrimerF                 PrimerR                 elevation   transect_name   ph
control     ATGACCAGATTA    GTGCCAGCMGCCGCGGTAA     GACTACHVHHHTWTCTAAT     control     control         control
BAQ4166.1   CAACTAGACTCG    GTGCCAGCMGCCGCGGTAA     GACTACHVHHHTWTCTAAT     4166        Baquedano       7.22
BAQ4166.2   GGAACGACGTGA    GTGCCAGCMGCCGCGGTAA     GACTACHVHHHTWTCTAAT     4166        Baquedano       7.22
BAQ4166.3   TGTCAGCTGTCG    GTGCCAGCMGCCGCGGTAA     GACTACHVHHHTWTCTAAT     4166        Baquedano       7.22
BAQ4697.1   CTGGTGCTGAAT    GTGCCAGCMGCCGCGGTAA     GACTACHVHHHTWTCTAAT     4697        Baquedano       7.44
BAQ4697.2   GACAGAGGTGCA    GTGCCAGCMGCCGCGGTAA     GACTACHVHHHTWTCTAAT     4697        Baquedano       7.44
BAQ4697.3   TCAGACCAACTG    GTGCCAGCMGCCGCGGTAA     GACTACHVHHHTWTCTAAT     4697        Baquedano       7.44
YUN3428.1   TACGCCCATCAG    GTGCCAGCMGCCGCGGTAA     GACTACHVHHHTWTCTAAT     3428        Yungay          7.20
YUN3428.2   AAGATCGTACTG    GTGCCAGCMGCCGCGGTAA     GACTACHVHHHTWTCTAAT     3428        Yungay          7.20
```

To avoid errors during the analysis, be careful to give the same samples between manifest and metadata file, and as possible, in the same order.

## Process parameters
This custom.config file control each process and it's parameters, **so it's very important to fulfill this file very carefully**. Otherwise, it's can lead to computational errors or misinterpretation of biological results.

In this section, we will describe the most important parameters for each process.

#### Main
```projectName```: the name of your project, **without space, tabulation or accented characters**

```input_manifest```: the full path to your manifest file. 

- ```input_manifest = "/path/to/samba/analysis/q2_manifest"```

```input_metadata```: the full path to your metadata file.
- ```input_metadata = "/path/to/samba/analysis/q2_metadata"```

#### Cleaning primers step using cutadapt
```primerF```: must be entered in the 5'-to-3' direction

```primerR```: must be entered in the 5'-to-3' direction

```overlap```: length shortest primer - 1

#### ASV taxonomic assignation using QIIME2 RDP-like program
```database```:

#### Decontamination step using microDecon package
```control_list```:

```nb_controls```:

```nb_samples```:

#### Taxonomy filtering
```tax_to_exclude```:

#### Differential abundance testing with ANCOM
```ancom_var```: 

#### Functional predictions with PICRUSt2
```picrust_var```: 

#### Remove samples for statistical analyses
```sample_to_remove```:

#### Statistics steps parameters
```alpha_div_group```: 

```beta_div_var```:
 
```desc_comp_crit```: 
