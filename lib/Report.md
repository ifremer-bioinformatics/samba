# <center><span style='text-align:center ; color:blue'>Report of your data processing using the QIIME 2 snakemake workflow created by the SeBiMER</span></center>



<span style="color:black"><u>**Authors:**</u> Laure QUINTRIC *(laure.quintric@ifremer.fr)* & Cyril NOËL *(cyril.noel@ifremer.fr)*  - SeBiMER (Ifremer)</span> 

<p align="center"> <img src="https://w3z.ifremer.fr/var/storage/images/_aliases/logo_main/medias-ifremer/banques-images-docs-bioinfo/cellule_bioinfo_2/1428369-1-fre-FR/Cellule_bioinfo_2.png" width="150"</p><br>


**Report established on Septembre 10, 2019 at 11:08**

[TOC]

## <span style="color:red">I. Report introduction</span>

<div style = "text-align: justify">
   This file was created using <b>Typora</b>, a free markdown editor and reader 
   (<i>https://typora.io/</i>). It allows to compile, in the form of a report, all results
   obtained by the metabarcoding analysis of your project 
   (<font color="blue">Risk Manche</font>) performed using <b>the QIIME 2 snakemake 
   workflow</b> (<i>qiime2_snakemake.sh</i>) created by the engineers of the 
   Ifremer's Bioinformatics department (<b>SeBiMER</b>).<br><br>
   <u>NB:</u> This report includes links. In order not to break them, be sure to have 
   extracted the entire compressed file that was sent to you and open this html file only
   in this folder. If you need to transfer the report, you should transfer the entire 
   directory.
</div>



## <span style="color:red">II. Description of the workflow</span>

<div style = "text-align: justify">
   The <i>qiime2_snakemake.sh</i> workflow developed by <u>Laure Quintric and 
   Cyril Noël</u> allows you to <b>computerize your metabarcoding analysis </b> using the 
   <b>DADA2 package through QIIME 2</b> (Bolyen <i>et al.</i>, 2019 ; version 2019.07). It
   also performs <b>basic diversity analyses</b> (alpha and beta) via the use of a 
   <font color="green">homemade R script</font> using mainly the vegan, phyloseq and 
   ggplot2 R packages. Specific statistics based on your metadata can also be carried out.
   This script is based on the use of the <b>Snakemake</b> workflow manager. All analysis 
   parameters are fully configurable via the <font color="red">config.yml file</font> 
   available in the workflow folder.
</div>



## <span style="color:red">III. Bioinformatic process</span>

### III.1. Importing sequencing data

<div style = "text-align: justify">
   The first step was to import sequencing data <b>under a QIIME 2 specific format</b>. 
   This was achieved by this step using the commands below. In addition, <b>descriptive 
   statistics of your data have been generated</b>.
</div> 

```bash
## Command run by snakemake :
### q2import.sh {manifest} {output.process} {output.visu} {output.summary}

# Import
qiime tools import \
--input-path {manifest} \
--output-path {output.process} \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-format PairedEndFastqManifestPhred33V2

# .qza to .qzv
qiime demux summarize --verbose \
--i-data  {output.process} \
--o-visualization {output.visu}

# Descriptive statistics
qiime tools export \
--input-path {output.visu} \
--output-path {output.summary}
```



<div style="background-color:rgba(135, 206, 250, 0.6)"><b><center>OUTPUT</center></b></div><br>
<ul style = "margin: 0 ; padding: 1 ; text-align: justify">
<li><span style="color:green"><b>step_01_q2import/data.qza</b> :</span> contains all data
   </li>
<li><span style="color:green"><b>step_01_q2import/data.qzv</b> :</span> contains all data 
   in a viewable format using QIIME 2 View (<i>https://view.qiime2.org/</i>)
   </li>
<li><span style="color:green"><b>step_01_q2import/summary_import</b> <i>[directory]</i> :
   </span> contains all descriptive statistics of your data
   </li>
</ul>



<div style="background-color:yellow";><b><center>RESULTS</center></b></div><br>
<ul style = "margin: 0 ; padding: 1 ; text-align: justify">
   <li>A total of 1,788,208 reads were imported</li>
   <li>Demultiplexed sequence counts summary :</li>
</ul>

|                    | Minimum | Median |  Mean  | Maximum |
| :----------------- | :-----: | :----: | :----: | :-----: |
| **Sequence count** |  4,816  | 32,168 | 33,115 | 59,441  |

<p align="center"> <img src="summary_import/demultiplex-summary.png" width="600"</p><br>
<ul style = "margin: 0 ; padding: 1 ; text-align: justify">
   <li>All descriptive statistics of your samples are available <a href="./summary_import/overview.html"><div style="display:inline-block;color:blue;"">here</div></a> (html output)</li>
</ul>



<div style="background-color:orange";><b><center>CONCLUSIONS</center></b></div><br>
<ul style="margin: 0 ; padding: 0 ; text-align: justify">
   No error was detected at this step. The rest of the workflow has been continued without
   any modification
</ul>




### III.2. Remove primers

<div style = "text-align: justify">
   The second step was to remove primers from sequences using <b>the cutadapt plugin</b>
   available in QIIME 2.
</div><br>

<ul style = "border: 2px dashed #FF0000 ;margin:0;padding:1;"><b><span style ="color:red;">CONFIGURATION SETTINGS</b></span>
   <li>primerF : CCTACGGGNGGCWGCAG</li>
   <li>primerR : GACTACHVGGGTATCTAAKCC</li>
   <li>errorRate : 0.1</li>
   <li>overlap : 16</li>
</ul><br>

```bash
## Command run by snakemake :
### q2_cutadapt.sh {resources.ncpus} {input} {params.primerF} {params.primerR} {params.errorRate} {params.overlap} {output.process} {output.visu} {output.summary}

# Trimming of the primers (run cutadapt)
qiime cutadapt trim-paired --verbose \
--p-cores {resources.ncpus} \
--i-demultiplexed-sequences {input} \
--p-front-f {params.primerF} \
--p-front-r {params.primerR} \
--p-error-rate {params.errorRate} \
--p-discard-untrimmed \
--p-match-read-wildcards \
--p-overlap {params.overlap} \
--o-trimmed-sequences {output.process}

# Summarize counts per sample for all samples
# qzv output
qiime demux summarize --verbose \
--i-data {output.process} \
--o-visualization {output.visu}

# html output
qiime tools export \
--input-path {output.visu} \
--output-path {output.summary}
```



<div style="background-color:rgba(135, 206, 250, 0.6)"><b><center>OUTPUT</center></b></div><br>
<ul style = "margin: 0 ; padding: 1 ; text-align: justify">
<li><span style="color:green"><b>step_02_q2cutadapt/data_trimmed.qza</b> :</span> contains
   all data whose sequence of the primers has been found and trimmed
   </li>
<li><span style="color:green"><b>step_02_q2cutadapt/data_trimmed.qzv</b> :</span> contains
   all data in a viewable format using QIIME 2 View (<i>https://view.qiime2.org/</i>)
   </li>
<li><span style="color:green"><b>step_02_q2cutadapt/summary_cutadapt</b> <i>[directory]</i>
   :</span> contains all descriptive statistics of your data
   </li>
</ul>



<div style="background-color:yellow";><b><center>RESULTS</center></b></div><br> 
<ul style = "margin: 0 ; padding: 1 ; text-align: justify">
   <li>A total of 891,344 reads were conserved after primer trimming. At this stage, 
       50.15% of reads eliminated 
   </li>
   <li>Demultiplexed sequence counts summary :</li>
</ul>

|                    | Minimum | Median |  Mean  | Maximum |
| :----------------- | :-----: | :----: | :----: | :-----: |
| **Sequence count** |  2,411  | 16,009 | 16,506 | 28,896  |

<p align="center"> <img src="summary_cutadapt/demultiplex-summary.png" width="600"</p><br>
<ul style = "margin: 0 ; padding: 1 ; text-align: justify">
   <li>All descriptive statistics of your samples after the trimming are available 
       <a href="./summary_cutadapt/overview.html"><div style="display:inline-block;
       color:blue;">here</div></a> (html output)
       </li>
   <li>The quality of your data before any quality filtering step is viewable
       <a href="./summary_cutadapt/quality-plot.html"><div style="display:inline-block;
       color:blue;">here</div></a> (html output)
   </li>
</ul>



<div style="background-color:orange";><b><center>CONCLUSIONS</center></b></div><br>
<ul style = "margin: 0 ; padding: 0 ; text-align: justify">
   No error was detected at this step. But, <u>a significant number of sequences have been
   eliminated</u>. Looking in detail, it turns out that: <font color="red"><b>1) forward 
   files are contaminated by reverse sequences</b></font> and <font color="red"><b>2) all
   files appear to be from two different sequencing</b></font>. In addition, the results 
   highlighted the duplication of barcode for different samples.<br>
   However, no solution for the moment has been identified. So, the rest of the workflow
   has still been continued without any modification
</ul>



### III.3.  Sequence quality control and feature table construction

<div style="text-align: justify">
   The third step of the workflow consisted in filtering the quality of the sequences 
   according to the <i>trim3F</i>, <i>trim3R</i>, <i>trunclenF</i>, <i>trunclenR</i>, 
   <i>minqual</i> and <i>maxee</i> parameters present in the config file. 
   Also to assemble the forward and reverse sequences, and to identify and remove the 
   chimeras according to the method defined by the <i>chimeras</i> parameter of the config
   file. To performed this step, the <b>DADA2 R package</b> was used through QIIME 2. 
</div><br>

<ul style = "border: 2px dashed #FF0000 ;margin:0;padding:1;"><b><span style ="color:red;">CONFIGURATION SETTINGS</b></span>
   <li>trim3F : 0</li>
   <li>trim3R : 0</li>
   <li>trunclenF : 260</li>
   <li>trunclenR : 220</li>
   <li>maxee : 2</li>
   <li>minqual : 2</li>
   <li>chimeras : consensus</li>
</ul>

```bash
## Command run by snakemake :
### q2_dada2.sh {input.data} {params.trim3F} {params.trim3R} {params.trunclenF} {params.trunclenR} {params.maxee} {params.minqual} {params.chimeras} {resources.ncpus} {output.process_repseqs} {output.process_table} {output.process_stats} {output.visu_stats} {output.visu_table} {input.metadata} {output.visu_repseqs} {output.export} {output.biom}

#Run dada2 : denoises paired-end sequences, dereplicates them and filters chimeras
qiime dada2 denoise-paired --verbose \
--i-demultiplexed-seqs {input.data} \
--p-trim-left-f {params.trim3F} --p-trim-left-r {params.trim3R} \
--p-trunc-len-f {params.trunclenF} --p-trunc-len-r {params.trunclenR} \
--p-max-ee {params.maxee} --p-trunc-q {params.minqual} \
--p-chimera-method {params.chimeras} --p-n-threads  {resources.ncpus} \
--o-representative-sequences {output.process_repseqs} \
--o-table {output.process_table} \
--o-denoising-stats {output.process_stats}

#Generate qzv file of the descriptive statistics of the dada2 process and extract results
qiime metadata tabulate --verbose \
--m-input-file {output.process_stats} \
--o-visualization {output.visu_stats}

qiime tools export \
--input-path {output.visu_stats} \
--output-path {output.export}

mv {output.export}/index.html {output.export}/STATS_DADA2.html 
mv {output.export}/metadata.tsv {output.export}/STATS_DADA2.tsv

#Generate qzv file of the feature table and extract results
qiime feature-table summarize --verbose \
--i-table {output.process_table} \
--o-visualization {output.visu_table} \
--m-sample-metadata-file {input.metadata}

qiime tools export \
--input-path {output.visu_table} \
--output-path {output.export}

#Generate qzv file of the representative sequences and extract them
qiime feature-table tabulate-seqs --verbose \
--i-data {output.process_repseqs} \
--o-visualization {output.visu_repseqs}

qiime tools export \
--input-path {output.visu_repseqs} \
--output-path {output.export}

mv {output.export}/sequences.fasta {output.export}/REPRESENTATIVE_SEQUENCES.fasta

#Generate the biom file of the ASV abundance (without taxonomy)
qiime tools export \
--input-path {output.process_table} \
--output-path {output.biom}
```



<div style="background-color:rgba(135, 206, 250, 0.6)"><b><center>OUTPUT</center></b></div><br>
<ul style = "margin: 0 ; padding: 1 ;text-align: justify">
   <li><font color="green"><b>step_03_q2dada2/rep_seqs.qza (.qzv)</b> :</font> contain 
       the representative sequences of the ASV (+viewable file)
   </li>
   <li><font color="green"><b>step_03_q2dada2/export_dada2/REPRESENTATIVE_SEQUENCES.fasta</b> :</font> 
                                 fasta file containing the representative sequence of the ASVs
   </li>
   <li><font color="green"><b>step_03_q2dada2/table.qza (.qzv)</b> :</font> contain 
       the abundance table (+viewable file)
   </li>
   <li><font color="green"><b>step_03_q2dada2/feature_table/feature-table.biom</b> :</font>
       the ASV abundance table in biom format (used after in the workflow)
   </li>
   <li><font color="green"><b>step_03_q2dada2/stats.qza (.qzv)</b> <i>[directory]</i> :</font> 
       contain the descriptive summary of your data obtained after this step with in 
       particular the number of sequences eliminated at each substep of Dada2 
       (quality filtering, denoising, merging and chimeras removal)
   </li>
   <li><font color="green"><b>step_03_q2dada2/export_dada2/STATS_DADA2.tsv</b> :</font> 
       tabulated file of the descriptive summary
   </li>
</ul>



<div style="background-color:yellow";><b><center>RESULTS</center></b></div><br>
<ul style = "margin: 0 ; padding: 1 ;text-align: justify">
   <li>The dynamics of the different step of filtering can be visualized in this 
       <a href="./export_dada2/STATS_DADA2.html"><div style="display:inline-block;
       color:blue;">html file</div></a> and are also available in a 
       <a href="./export_dada2/STATS_DADA2.tsv"><div style="display:inline-block;
       color:blue;">.tsv</div></a>
   </li>
</ul>

| DADA2 step  |  Input  |    Filtered     |    Denoised     |     Merged      |  Non-chimeric   |
| :---------: | :-----: | :-------------: | :-------------: | :-------------: | :-------------: |
| Reads count | 891,344 | 510,865 (42.7%) | 407,370 (54.3%) | 255,952 (71.3%) | 230,955 (74.1%) |
<ul style = "margin: 0 ; padding: 1 ;text-align: justify">
   <li>Distribution of sequences in samples
   </li>
</ul>  
<p align="center"> <img src="export_dada2/sample-frequencies.png" width="600"</p><br>

| Min ASV count | Mean ASV count | Median ASV count | Max ASV count |
| :-----------: | :------------: | :--------------: | :-----------: |
|      318      |     4,277      |      4,133       |     8,985     |

<ul style = "margin: 0 ; padding: 1 ;text-align: justify">
   <li>Feature details are available <a href="./export_dada2/feature-frequency-detail.html">
       <div style="display:inline-block;color:blue;">here</div></a> (html output)
   </li>
</ul>

| Sequence count | Min Length | Max Length | Mean Length | Range | Standard Deviation |
| :------------: | :--------: | :--------: | :---------: | :---: | :----------------: |
|     2,476      |    318     |    436     |   419.13    |  118  |       10.39        |

| Percentile:  |  2%  |  9%  | 25%  | 50%  | 75%  | 91%  | 98%  |
| :----------: | :--: | :--: | :--: | :--: | :--: | :--: | :--: |
| Length (nts) | 402  | 402  | 407  | 422  | 427  | 427  | 428  |

<ul style = "margin: 0 ; padding: 1 ;text-align: justify">
   <li>Details about the samples can be found by going to 
       <a href="./export_dada2/sample-frequency-detail.html">
       <div style="display:inline-block;color:blue;">this interactive html page</div></a>
   </li>
   <li>Finally, you can retrieved the reference sequences of your ASVs in 
       <a href="./export_dada2/REPRESENTATIVE_SEQUENCES.fasta">
       <div style="display:inline-block;color:blue;">this fasta</div></a>
   </li>
</ul> 



<div style="background-color:orange";><b><center>CONCLUSIONS</center></b></div><br>
<ul style = "margin: 0 ; padding: 0 ;text-align: justify">
   No error was detected at this step. But, <u>a significant number of sequences have been 
   eliminated</u> mainly due to their <font color="red"><b>low quality</b></font> that did not fit into the parameters used 
   or after the trimming could not be joined.
</ul>



### III.4.  Taxonomic assignation

<div style="text-align: justify">
   The taxonomic assignation performed during the fourth step of the workflow allowed to 
   affiliate each ASV to a taxonomy by using as reference the database defined in the 
   config file 
</div><br>

<ul style = "border: 2px dashed #FF0000 ;margin:0;padding:1;"><b><span style ="color:red;">CONFIGURATION SETTINGS</b></span>
   <li>database : /home/ref-bioinfo/tools_data/qiime2/DATABASE_silva_99_16S.qza</li>
   <li>confidence : 0.8</li>
</ul>



```bash
## Command run by snakemake :
### q2_taxonomy.sh {resources.ncpus} {params.confidence} {params.database} {input.repseqs} {output.taxonomy} {output.visu} {output.export} {output.tax_tsv}

#Taxonomy assignment
qiime feature-classifier classify-sklearn \
--p-n-jobs {resources.ncpus} \
--p-confidence {params.confidence} \
--i-classifier {params.database} \
--i-reads {input.repseqs} \
--o-classification {output.taxonomy}

#Convert the taxonomy into a qzv file (visualizable with Qiime2view) and extract results
qiime metadata tabulate \
--m-input-file {output.taxonomy} \
--o-visualization {output.visu}

qiime tools export \
--input-path {output.visu} \
--output-path {output.export}

#Rename tabular taxonomy file and modify header
mv {output.export}/metadata.tsv {output.tax_tsv}
sed -i '1,2d' {output.tax_tsv}
sed -i '1 i\#OTUID	taxonomy	confidence' {output.tax_tsv}
```



<div style="background-color:rgba(135, 206, 250, 0.6)"><b><center>OUTPUT</center></b></div><br>
<ul style = "margin: 0 ; padding: 1 ;text-align: justify">
   <li><font color="green"><b>step_04_q2taxonomy/taxonomy.qza (.qzv)</b> :</font> contain 
       the taxonomy affiliated at each ASV
   </li>
   <li><font color="green"><b>step_04_q2taxonomy/ASV_taxonomy.tsv</b> :</font> tabulated 
       file containing the taxonomy and formatted for subsequent use
   </li>
</ul>



<div style="background-color:yellow";><b><center>RESULTS</center></b></div><br>
<ul style = "margin: 0 ; padding: 1 ;text-align: justify">
   <li>The result of the taxonomic affiliation is available in 
       <a href="./export_taxonomy/index.html"><div style="display:inline-block;color:blue;">
       html</div></a> and <a href="./ASV_taxonomy.tsv"><div style="display:inline-block;
       color:blue;">tabulated</div></a> formats
   </li>
</ul> 



<div style="background-color:orange";><b><center>CONCLUSIONS</center></b></div><br>
<ul style="margin: 0 ; padding: 0 ;text-align: justify">
   No error was detected at this step
</ul>



### III.5.  Final output

<div style="text-align: justify">
   This is the last step of the bioinformatic process where the goal was to merge the ASV 
   abundance table with the taxonomy file
</div><br/>

```bash
## Command run by snakemake :
### q2_output.sh {input.table} {input.taxonomy} {output.biom} {output.tsv}

#Add taxonomy to biom ASV abundance table
biom add-metadata -i biom add-metadata \
-i {input.table} \
--observation-metadata-fp {input.taxonomy} \
-o {output.biom} \
--sc-separated taxonomy

# Convert into a tabulated format
biom convert \
-i {output.biom} \
-o {output.tsv} \
--to-tsv --header-key taxonomy
```



<div style="background-color:rgba(135, 206, 250, 0.6)"><b><center>OUTPUT</center></b></div><br>
<ul style = "margin: 0 ; padding: 1 ;text-align: justify">
   <li><font color="green"><b>step_05_q2output/Final_ASV_table_with_taxonomy.biom (.tsv)</b> :</font> 
       final output in biom and tabulated format
   </li>
</ul>



<div style="background-color:yellow";><b><center>RESULTS</center></b></div><br>
<ul style = "margin: 0 ; padding: 1 ;text-align: justify">
   <li>The final output is available in <a href="./Final_ASV_table_with_taxonomy.biom">
       <div style="display:inline-block;color:blue;">biom</div></a> and 
       <a href="./Final_ASV_table_with_taxonomy.tsv"><div style="display:inline-block;
       color:blue;">tabulated</div></a> formats
   </li>
</ul> 



<div style="background-color:orange";><b><center>CONCLUSIONS</center></b></div><br>
<ul style="margin: 0 ; padding: 0 ;text-align: justify">
   No error was detected at this step
</ul>



## <span style="color:red">IV. General statistical analyses</span>

<ul style = "border: 2px dashed #FF0000 ;margin:0;padding:1;"><b><span style ="color:red;">CONFIGURATION SETTINGS</b></span>
   <li>project_name : RiskManche_cDNA</li>
   <li>perc_abund_threshold : 1</li>
   <li>distance : bray</li>
   <li>column_sample_replicat : Sample_Location</li>
</ul>

### IV.1.  Alpha diversity

<p align="center"> <img src="./R/FIGURES/alpha_diversity/alpha_div_plots.svg"></p><br>
### IV.2.  Taxonomic barplots

<p align="center"> <img src="./R/FIGURES/alpha_diversity/barplot_relabund_phylum.svg"></p><br>
<p align="center"> <img src="./R/FIGURES/alpha_diversity/barplot_relabund_family.svg"></p><br>
<p align="center"> <img src="./R/FIGURES/alpha_diversity/barplot_relabund_genus.svg"></p><br>
### IV.3. Heatmap representation

<p align="center"> <img src="./R/FIGURES/alpha_diversity/heatmap_class.svg"></p><br>
<p align="center"> <img src="./R/FIGURES/alpha_diversity/heatmap_family.svg"></p><br>
<p align="center"> <img src="./R/FIGURES/alpha_diversity/heatmap_genus.svg"></p><br>
### IV.4.  Beta diversity

#### <span style = "color:grey;"><center>IV.4.a. Non-normalized data</center></span>

##### <i>IV.4.a.i. ASV resolution</i>

<p align="center"> <img src="./R/FIGURES/beta_diversity/ASV_ordination_plot.svg"></p><br>
<p align="center"> <img src="./R/FIGURES/beta_diversity/ASV_ordination_plot_wrapped.svg"></p><br>
##### <i>IV.4.a.ii. Sample resolution</i>

<p align="center"> <img src="./R/FIGURES/beta_diversity/samples_ordination_plot.svg"></p><br>
##### <i>IV.4.a.iii. Split graphic</i>

<p align="center"> <img src="./R/FIGURES/beta_diversity/split-graph_ordination_plot.svg"></p><br>
#### <span style = "color:grey;"><center>IV.4.b. Rarefied data</center></span>

##### <i>IV.4.b.i. ASV resolution</i>

<p align="center"> <img src="./R/FIGURES/beta_diversity_rarefied/ASV_ordination_plot_rarefied.svg"></p><br>
<p align="center"> <img src="./R/FIGURES/beta_diversity_rarefied/ASV_ordination_plot_wrapped_rarefied.svg"></p><br>
##### <i>IV.4.b.ii. Sample resolution</i>

<p align="center"> <img src="./R/FIGURES/beta_diversity_rarefied/samples_ordination_plot_rarefied.svg"></p><br>
##### <i>IV.4.b.iii. Split graphic</i>

<p align="center"> <img src="./R/FIGURES/beta_diversity_rarefied/split-graph_ordination_plot_rarefied.svg"></p><br>
#### <span style = "color:grey;"><center>IV.4.c. Data normalized using DESeq2</center></span>

##### <i>IV.4.c.i. ASV resolution</i>

<p align="center"> <img src="./R/FIGURES/beta_diversity_DESeq2/ASV_ordination_plot_DESeq2.svg"></p><br>
<p align="center"> <img src="./R/FIGURES/beta_diversity_DESeq2/ASV_ordination_plot_wrapped_DESeq2.svg"></p><br>
##### <i>IV.4.c.ii. Sample resolution</i>

<p align="center"> <img src="./R/FIGURES/beta_diversity_DESeq2/samples_ordination_plot_DESeq2.svg"></p><br>
##### <i>IV.4.c.iii. Split graphic</i>

<p align="center"> <img src="./R/FIGURES/beta_diversity_DESeq2/split-graph_ordination_plot_DESeq2.svg"></p><br>
#### <span style = "color:grey;"><center>IV.4.d. Data normalized using CSS (metagenomeSeq R package)</center></span>

##### <i>IV.4.d.i. ASV resolution</i>

<p align="center"> <img src="./R/FIGURES/beta_diversity_CSS/ASV_ordination_plot_CSS.svg"></p><br>
<p align="center"> <img src="./R/FIGURES/beta_diversity_CSS/ASV_ordination_plot_wrapped_CSS.svg"></p><br>
##### <i>IV.4.d.ii. Sample resolution</i>

<p align="center"> <img src="./R/FIGURES/beta_diversity_CSS/samples_ordination_plot_CSS.svg"></p><br>
##### <i>IV.4.d.iii. Split graphic</i>

<p align="center"> <img src="./R/FIGURES/beta_diversity_CSS/split-graph_ordination_plot_CSS.svg"></p><br>
## <span style="color:red">V. Statistical analyses specific to the project</span>

### 
