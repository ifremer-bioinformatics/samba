# Training Dataset
## Training dataset used from [Qiime2 Tutorial] (https://docs.qiime2.org/2019.7/tutorials/atacama-soils/), , [associated publication](https://msystems.asm.org/content/2/3/e00195-16)

In this directory, you can find : 

* a directory named **dna-sequence-raw** containing: 
    * forward (R1) and reverse (R2) fastq files from 46 samples 

* a **q2_manifest** file required by QIIME 2 and containing:
    * 3 columns representing for each sample -> the name of the sample, the path to the R1 file and the path to the R2 file
    
* a **q2_metadata** file required by QIIME 2 and for statistical analyses. It contains:
    * a first column named *sampleid* -> name of the sample (must be identical to the sampleid in the q2_manifest)
    * a second column named *barcode* -> barcode of each sample
    * the third and fourth columns represent the forward and reverse primers used in the study (V4 region of the 16S)
    * others columns represent variables describing each sample used for the statistical analyses  

## How to use this test dataset

Simply run the script [RunSAMBA_training_dataset.sh](https://gitlab.ifremer.fr/bioinfo/SAMBA-nextflow/blob/master/RunSAMBA_training_dataset.sh) directly from the root of SAMBA

```
./RunSAMBA_training_dataset.sh
```
