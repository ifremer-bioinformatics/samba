# <center><span style='text-align:center ; color:blue'>Report of your SAMBA metabarcoding analysis</span></center>



<span style="color:black"><u>**Authors:**</u> Laure QUINTRIC & Cyril NOËL - Bioinformatics engineers  - SeBiMER (Ifremer)</span> 

<span style="color:black"><u>**Contact:**</u> *samba-sebimer@ifremer.fr*</span> 

[TOC]

## <span style="color:red">I. Report introduction</span>

<div style = "text-align: justify">
   This document is the report of your metabarcoding analysis performed using using 
   <b>the SAMBA workflow</b> created by the engineers of the Ifremer's Bioinformatics 
   service (<b>SeBiMER</b>).
   This report includes results in the form of clickable links referring to generated 
   analysis files. In order not to break them, <b>be sure not to move this document 
   independently of all the result folder</b>.
</div>



## <span style="color:red">II. Description of the workflow</span>

<div style = "text-align: justify">
   The <i>SAMBA workflow</i> developed by <u>Laure Quintric and Cyril Noël</u> allows 
   you to <b>computerize metabarcoding analysis </b> (16S, 18S, ITS) using the 
   <b>DADA2 package through QIIME 2</b> (Bolyen <i>et al.</i>, 2019 ; version 2019.10). 
   It also performs <b>basic diversity analyses</b> (alpha and beta) via the use of 
   <font color="green">homemade R scripts</font> using the vegan, phyloseq and ggplot2 
   R packages. 
   Specific statistics based on your metadata can also be carried out either by using, 
   in your own scripts, the <i>phyloseq</i> object created by SAMBA (advanced users) or 
   by coming back to us for additional analyses.
   This workflow is based on the use of the <b>NextFlow</b> workflow managers. All 
   analysis parameters are fully configurable via the 
   <font color="red">params.config file</font> available in the config folder of SAMBA.
</div>



## <span style="color:red">III. Bioinformatic process</span>

### III.1. Data integrity

<div style="text-align: justify"><span style="color:black">
   This first step allows to <b>analyze the integrity of your data</b> in order to identify 
   potential problems related to the sequencing processes. It checks:
<ul>
    <li>that each read is correctly associated with the proper sample (sequence barcode verification)</li>
    <li>that they come from a single sequencer</li>
    <li>the efficiency of forward and reverse PCR amplification</li>
</ul>
   <b>This step is carried out by a <font color="green">homemade script</font></b>
</span></div>



<div style="background-color:yellow"><center><b>RESULTS</b></center></div><br>
<ul style = "margin: 0 ; padding: 1 ;">
   <li>All data integrity verification are summarized in <a href="./data_integrity.csv"><div style="display:inline-block;color:blue;"">this csv file</div></a></li>
</ul>



### III.2. Importing raw data

<div style="text-align: justify"><span style="color:black">
   The step performs the import of sequencing data directly from DATAREF <b>into a QIIME 2 
   specific format</b>. In addition, <b>descriptive statistics of your data are generated
   </b>. The SAMBA workflow ran for you the commands below :
</span></div>

```bash
### Import data using the q2_manifest ###

qiime tools import \
--input-path q2_manifest \
--output-path data.qza \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-format PairedEndFastqManifestPhred33V2

### Descriptive statistics ###

# Summarize statistics
qiime demux summarize \
--verbose \
--i-data data.qza \
--o-visualization data.qzv

# Export statistics
qiime tools export \
--input-path data.qzv \
--output-path import_output
```



<div style="background-color:rgba(135, 206, 250, 0.6)"><b><center>OUTPUT FOLDER: <a href="./import_output"><div style="display:inline-block;color:blue;"">import_output</div></a></center></b>
</div>



<div style="background-color:yellow"><center><b>RESULTS</b></center></div><br>
<ul style = "margin: 0 ; padding: 1 ;">
   <li>Sample repartition according to the sequence count</li>
</ul>  
<p align="center"> <img src="./import_output/demultiplex-summary.png" width="600"</p><br>
<ul style = "margin: 0 ; padding: 1 ;">
   <li>All descriptive statistics of your samples are available <a href="./import_output/overview.html"><div style="display:inline-block;color:blue;"">here</div></a> (html output)</li>
</ul>




### III.3. Removal of primers

<div style="text-align: justify"><span style="color:black">
   The third step is to remove the primers using <b>the cutadapt plugin</b> available in 
   QIIME 2. In addition, <b>descriptive statistics of your data are generated</b>. 
   The following commands made this step possible using the parameters that you have 
   defined in the parameters configuration file :
</span></div><br>

<ul style = "border: 2px dashed #FF0000 ;margin:0;padding:1;"><b><span style ="color:red;">CONFIGURATION SETTINGS USED</b></span>
   <li>primerF</li>
   <li>primerR</li>
   <li>errorRate</li>
   <li>overlap</li>
</ul><br>

```bash
### Cutadapt process ###

qiime cutadapt trim-paired \
--verbose \
--p-cores 14 \
--i-demultiplexed-sequences data.qza \
--p-front-f $primerF \
--p-front-r $primerR \
--p-error-rate $errorRate \
--p-discard-untrimmed \
--p-match-read-wildcards \
--p-overlap $overlap \
--o-trimmed-sequences data_trimmed.qza

### Descriptive statistics ###

# Summarize statistics of the cutadapt process
qiime demux summarize \
--verbose \
--i-data data_trimmed.qza \
--o-visualization data_trimmed.qzv

# Export statistics
qiime tools export \
--input-path data_trimmed.qzv \
--output-path trimmed_output   
```



<div style="background-color:rgba(135, 206, 250, 0.6)"><b><center>OUTPUT FOLDER: <a href="./trimmed_output"><div style="display:inline-block;color:blue;"">trimmed_output</div></a></center></b>
</div>



<div style="background-color:yellow";><b><center>RESULTS</center></b></div><br> 
<ul style = "margin: 0 ; padding: 1 ;">
   <li>Sample repartition according to the sequence count</li></ul>
<p align="center"> <img src="./trimmed_output/demultiplex-summary.png" width="600"</p><br>
<ul style = "margin:0 ; padding: 1 ;">
   <li>All descriptive statistics of your samples after the trimming are available <a href="./trimmed_output/overview.html"><div style="display:inline-block;color:blue;">here</div></a> (html output)</li>
   <li>The quality of your data before any quality filtering step is viewable <a href="./trimmed_output/quality-plot.html"><div style="display:inline-block;color:blue;">here</div></a> (html output)</li>
</ul>   



### III.4.  Sequence quality control and feature table construction

<div style="text-align: justify"><span style="color:black">
   The following step of the workflow is to filter the quality of the sequences according
   to the parameters defined yourself by having considered the quality of your data in the
   previous step. 
   Also to assemble the forward and reverse sequences, and to identify and remove the 
   chimeras. To performed this step, the <b>DADA2 R package</b> was used through QIIME 2 
   using the following commands :
</span></div><br>

<ul style = "border: 2px dashed #FF0000 ;margin:0;padding:1;"><b><span style ="color:red;">CONFIGURATION SETTINGS USED</b></span>
   <li>trim5F</li>
   <li>trim5R</li>
   <li>trunclenF</li>
   <li>trunclenR</li>
   <li>maxeeF</li>
   <li>maxeeR</li>
   <li>minqual</li>
   <li>chimeras</li>
</ul>



```bash
### All DADA2 processes ###
qiime dada2 denoise-paired \
--verbose \
--i-demultiplexed-seqs data_trimmed.qza \
--p-trim-left-f $trim5F \
--p-trim-left-r $trim5R \
--p-trunc-len-f $trunclenF \
--p-trunc-len-r $trunclenR \
--p-max-ee-f $maxeeF \
--p-max-ee-r $maxeeR \
--p-trunc-q $minqual \
--p-chimera-method $chimeras \
--p-n-threads 14 \
--o-representative-sequences rep_seqs.qza \
--o-table table.qza \
--o-denoising-stats stats.qza

### Descriptive statistics ###

# Summarize statistics of the DADA2 processes
qiime metadata tabulate \
--verbose \
--m-input-file stats.qza \
--o-visualization stats.qzv

# Export statistics
qiime tools export \
--input-path stats.qzv \
--output-path dada2_output

### ASV table ###

# Construction of the count table
qiime tools export \
--input-path table.qza \
--output-path dada2_output  

# Summarize statistics from the table file
qiime feature-table summarize \
--verbose \
--i-table table.qza \
--o-visualization table.qzv \
--m-sample-metadata-file q2_metadata

# Export of the table ASV statistics
qiime tools export \
--input-path table.qzv \
--output-path dada2_output

### Reference sequences ###

# Obtaining ASV reference sequences
qiime feature-table tabulate-seqs \
--verbose \
--i-data rep_seqs.qza \
--o-visualization rep_seqs.qzv

# Export of ASV reference sequences
qiime tools export \
--input-path rep_seqs.qzv \
--output-path dada2_output
```



<div style="background-color:rgba(135, 206, 250, 0.6)"><b><center>OUTPUT FOLDER: <a href="./dada2_output"><div style="display:inline-block;color:blue;"">dada2_output</div></a></center></b></div>



<div style="background-color:yellow";><b><center>RESULTS</center></b></div><br>
<ul style = "margin: 0 ; padding: 1 ;text-align: justify">
   <li>The dynamics of the different step of filtering can be visualized in this <a href="./dada2_output/index.html"><div style="display:inline-block;color:blue;">html file</div></a> and are also available in a <a href="./dada2_output/metadata.tsv"><div style="display:inline-block;color:blue;">tabulated file</div></a></li>
   <li>Distribution of sequences in samples</li>
        <p align="center"><img src="./dada2_output/sample-frequencies.png" width="600"</p><br>
   <li>Feature details are available <a href="./dada2_output/feature-frequency-detail.html"><div style="display:inline-block;color:blue;">here</div></a> (html output)</li>
   <li>Details about the samples can be found by going to <a href="./dada2_output/sample-frequency-detail.html"><div style="display:inline-block;color:blue;">this interactive html page</div></a></li>
           <li>Finally, you can retrieved the reference sequences of your ASVs in <a href="./dada2_output/sequences.fasta"><div style="display:inline-block;color:blue;">this fasta</div></a></li>
</ul>



### III.5.  Taxonomic assignation

<div style="text-align: justify">
   The taxonomic assignation performed during the fourth step of the workflow allowed to 
   affiliate each ASV to a taxonomy by using as reference the database defined by yourself
</div><br>

<ul style = "border: 2px dashed #FF0000 ;margin:0;padding:1;"><b><span style ="color:red;">CONFIGURATION SETTINGS USED</b></span>
   <li>database</li>
   <li>confidence</li>
</ul>



```bash
### Taxonomic assignation ###

qiime feature-classifier classify-sklearn \
--p-n-jobs 14 \
--p-confidence $confidence \
--i-classifier $database \
--i-reads rep_seqs.qza \
--o-classification taxonomy.qza

### Output formatting ###

# Summarize taxonomic results
qiime metadata tabulate \
--m-input-file taxonomy.qza \
--o-visualization taxonomy.qzv

# Export
qiime tools export \
--input-path taxonomy.qzv \
--output-path taxo_output

# Formatting
mv taxo_output/metadata.tsv ASV_taxonomy.tsv
sed -i '1,2d' ASV_taxonomy.tsv
sed -i '1 i\#OTUID\ttaxonomy\tconfidence' ASV_taxonomy.tsv
```



<div style="background-color:rgba(135, 206, 250, 0.6)"><b><center>OUTPUT FOLDER: <a href="./taxo_output"><div style="display:inline-block;color:blue;"">taxo_output</div></a></center></b></div>



<div style="background-color:yellow";><b><center>RESULTS</center></b></div><br>
<ul style = "margin: 0 ; padding: 1 ;text-align: justify">
   <li>The result of the taxonomic affiliation is available in <a href="./taxo_output/index.html"><div style="display:inline-block;color:blue;">html</div></a> and <a href="./taxo_output/ASV_taxonomy.tsv"><div style="display:inline-block;color:blue;">tabulated</div></a> formats</li>
</ul> 



### III.6.  Final outputs

<div style="text-align: justify"><span style="color:black">
   This is the last step of the bioinformatic process where the goal is to merge the ASV abundance table with the taxonomy file. This is performed using the following commands :
</span></div>

```bash
biom add-metadata \
-i dada2_output/feature-table.biom \
--observation-metadata-fp ASV_taxonomy.tsv \
-o Final_ASV_table_with_taxonomy.biom \
--sc-separated taxonomy

biom convert \
-i Final_ASV_table_with_taxonomy.biom \
-o Final_ASV_table_with_taxonomy.tsv \
--to-tsv \
--header-key taxonomy
```



<div style="background-color:yellow";><b><center>RESULTS</center></b></div><br>
<ul style = "margin: 0 ; padding: 1 ;text-align: justify">
   <li>The final ASV table can be viewed <a href="./Final_ASV_table_with_taxonomy.tsv"><div style="display:inline-block;color:blue;">here</div></a>. A biom file is available <a href="./Final_ASV_table_with_taxonomy.biom"><div style="display:inline-block;color:blue;">here</div></a> formats</li>
</ul> 



## <span style="color:red">IV. General statistical analyses</span>

<ul style = "border: 2px dashed #FF0000 ;margin:0;padding:1;"><b><span style ="color:red;">CONFIGURATION SETTINGS</b></span>
   <li>Normalization method : "None / Rarefaction / DESeq2 / CSS"</li>
   <li>Grouping variable</li>
   <li>Distance matrix : "Jaccard / Bray-Curtis / UniFrac / Weighted UniFrac"</li>
</ul>



### IV.1.  Alpha diversity

#### <i>IV.1.a diversity indices</i>

<p align="center"> <img src="./R/FIGURES/alpha_diversity/diversity_index/alpha_div_plots.png" width="1000"</p><br>

#### <i>IV.1.b taxonomic diversity</i>

<p align="center"> <img src="./R/FIGURES/alpha_diversity/diversity_barplots/barplot_phylum.png" width="1000"</p><br>
<p align="center"> <img src="./R/FIGURES/alpha_diversity/diversity_barplots/barplot_class.png" width="1000"</p><br>
<p align="center"> <img src="./R/FIGURES/alpha_diversity/diversity_barplots/barplot_order.png" width="1000"</p><br>
<p align="center"> <img src="./R/FIGURES/alpha_diversity/diversity_barplots/barplot_genus.png" width="1000"</p><br>
<ul style = "margin: 0 ; padding: 1 ;text-align: justify">
   <li>The results of the significance tests carried out for each variable on each diversity indexes can be viewed <a href="./R/FIGURES/alpha_diversity/index_significance_tests.txt"><div style="display:inline-block;color:blue;">here</div></a></li>
</ul> 



### IV.2.  Beta diversity

<div style="text-align: justify"><span style="color:black">
   <font color="red"><b>NOTE 1: Raw results without any data normalization are available 
      <a href="./R/FIGURES/beta_diversity_non_normalized"><div style="display:inline-block;color:blue;">here</div></a>. However, we do not recommend using them because from a statistical point of view, these results are very questionable.</b>
   </font><br><br>
   <font color="red"><b>NOTE 2: For the beta diversity, SAMBA calculates 4 different distance 
      matrices with jaccard, bray-curtis, UniFrac and weighted UniFrac indices. To facilitate 
      the reading of this report, we have displayed that the results of the bray-curtis and 
      weighted UniFrac matrices. All results are available in each normalization folder <a href="./R/FIGURES/"><div style="display:inline-block;color:blue;">here</div></a></b>
   </font>
</span></div>



#### <i>IV.2.a NMDS</i> 

###### IV.2.a.i rarefied data

<p align="center"> <img src="./R/FIGURES/beta_diversity_rarefied/NMDS/NMDS_rarefied_bray.png" width="1000"</p><br>
<p align="center"> <img src="./R/FIGURES/beta_diversity_rarefied/NMDS/NMDS_rarefied_wunifrac.png" width="1000"</p><br>

###### IV.2.a.ii DESeq2 normalization

<p align="center"> <img src="./R/FIGURES/beta_diversity_DESeq2/NMDS/NMDS_DESeq2_bray.png" width="1000"</p><br>
<p align="center"> <img src="./R/FIGURES/beta_diversity_DESeq2/NMDS/NMDS_DESeq2_wunifrac.png" width="1000"</p><br>

###### IV.2.a.iii CSS normalization

<p align="center"> <img src="./R/FIGURES/beta_diversity_CSS/NMDS/NMDS_CSS_bray.png" width="1000"</p><br>
<p align="center"> <img src="./R/FIGURES/beta_diversity_CSS/NMDS/NMDS_CSS_wunifrac.png" width="1000"</p><br>

#### <i>IV.2.b PCoA</i> 

###### IV.2.b.i rarefied data

<p align="center"> <img src="./R/FIGURES/beta_diversity_rarefied/PCoA/PCoA_rarefied_bray.png" width="1000"</p><br>
<p align="center"> <img src="./R/FIGURES/beta_diversity_rarefied/PCoA/PCoA_rarefied_wunifrac.png" width="1000"</p><br>

###### IV.2.b.ii DESeq2 normalization

<p align="center"> <img src="./R/FIGURES/beta_diversity_DESeq2/PCoA/PCoA_DESeq2_bray.png" width="1000"</p><br>
<p align="center"> <img src="./R/FIGURES/beta_diversity_DESeq2/PCoA/PCoA_DESeq2_wunifrac.png" width="1000"</p><br>

###### IV.2.b.iii CSS normalization

<p align="center"> <img src="./R/FIGURES/beta_diversity_CSS/PCoA/PCoA_CSS_bray.png" width="1000"</p><br>
<p align="center"> <img src="./R/FIGURES/beta_diversity_CSS/PCoA/PCoA_CSS_wunifrac.png" width="1000"</p><br>

#### <i>IV.2.c Hierarchical clustering</i> 

###### IV.2.c.i rarefied data

<p align="center"> <img src="./R/FIGURES/beta_diversity_rarefied/Hierarchical_Clustering/hclustering_rarefied_bray.png" width="1000"</p><br>
<p align="center"> <img src="./R/FIGURES/beta_diversity_rarefied/Hierarchical_Clustering/hclustering_rarefied_wunifrac.png" width="1000"</p><br>

###### IV.2.c.ii DESeq2 normalization

<p align="center"> <img src="./R/FIGURES/beta_diversity_DESeq2/Hierarchical_Clustering/hclustering_DESeq2_bray.png" width="1000"</p><br>
<p align="center"> <img src="./R/FIGURES/beta_diversity_DESeq2/Hierarchical_Clustering/hclustering_DESeq2_wunifrac.png" width="1000"</p><br>

###### IV.2.c.iii CSS normalization

<p align="center"> <img src="./R/FIGURES/beta_diversity_CSS/Hierarchical_Clustering/hclustering_CSS_bray.png" width="1000"</p><br>
<p align="center"> <img src="./R/FIGURES/beta_diversity_CSS/Hierarchical_Clustering/hclustering_CSS_wunifrac.png" width="1000"</p><br>

#### <i>IV.2.d Explained variance</i> 

###### IV.2.d.i rarefied data

<p align="center"> <img src="./R/FIGURES/beta_diversity_rarefied/pie_ExpVar_rarefied_bray.png" width="1000"</p><br>
<p align="center"> <img src="./R/FIGURES/beta_diversity_rarefied/pie_ExpVar_rarefied_wunifrac.png" width="1000"</p><br>

###### IV.2.d.ii DESeq2 normalization

<p align="center"> <img src="./R/FIGURES/beta_diversity_DESeq2/pie_ExpVar_DESeq2_bray.png" width="1000"</p><br>
<p align="center"> <img src="./R/FIGURES/beta_diversity_DESeq2/pie_ExpVar_DESeq2_wunifrac.png" width="1000"</p><br>

###### IV.2.d.iii CSS normalization

<p align="center"> <img src="./R/FIGURES/beta_diversity_CSS/pie_ExpVar_CSS_bray.png" width="1000"</p><br>
<p align="center"> <img src="./R/FIGURES/beta_diversity_CSS/pie_ExpVar_CSS_wunifrac.png" width="1000"</p><br>