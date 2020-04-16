#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Script name: data_integrity.sh                                          ####
##                                                                           ##
## Purpose of script: Check integrity of raw data                            ##
##                                                                           ##
##                                                                           ##
##                                                                           ##
## Authors: Laure QUINTRIC and Cyril NOEL                                  ####
##          Bioinformatics engineers                                         ##
##          SeBiMER, Ifremer                                                 ##
##                                                                           ##
## Creation Date: 2019-08-29                                               ####
## Modified on: 2020-01-31                                                 ####
##                                                                           ##
## Email: samba-sebimer@ifremer.fr                                         ####
##                                                                           ##
## Copyright (c) SeBiMER, august-2019                                      ####
## This program is free software: you can redistribute it and/or modify it   ##
## under the terms of the GNU Affero General Public License as published by  ##
## the Free Software Foundation, either version 3 of the License, or         ##
## (at your option) any later version.                                       ## 
##                                                                           ##
## License at https://www.gnu.org/licenses/agpl-3.0.txt                      ##
##                                                                           ##
## This program is distributed in the hope that it will be useful, but       ##
## WITHOUT ANY WARRANTY; without even the implied warranty of                ##
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                      ##
## See the GNU Affero General Public License for more details.               ##
##                                                                           ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Command run by nextflow :
## data_integrity.sh ${manifest} ${metadata} ${params.data_integrity.primerF} ${params.data_integrity.primerR} data_integrity.txt verifications.ok verifications.bad ${params.data_integrity.barcode} ${params.data_integrity.sampleid_column_name} ${params.data_integrity.R1_files_column_name} ${params.data_integrity.R2_files_column_name} ${params.data_integrity.barcode_filter} ${params.data_integrity.primer_filter} > data_integrity.log 2>&1

# Arguments
args=("$@")
manifest=${args[0]}
metadata=${args[1]}
primerF=${args[2]}
primerR=${args[3]}
summary=${args[4]}
verif_ok=${args[5]}
verif_bad=${args[6]}
barcode=${args[7]}
sampleid=${args[8]}

singleEnd=${args[14]}
if ${singleEnd}; then
    R1_files=${args[9]}
else 
    R1_files=${args[10]}
    R2_files=${args[11]}
fi

barcode_filter=${args[12]}
primer_filter=${args[13]}

# Verify if metadata file contains NA values
grep -P "NA\t|\tNA" ${metadata}
[ $? -eq 0 ] && { 
   echo "NA values found in $metadata, please remove them before running SAMBA"; 
   touch ${verif_bad} && echo "@@@ --- All verifications are not satisfied --- @@@"; 
   exit 0; 
}

# Creation of temporary files
#get sample id list from manifest file
COL=`head -n1 ${manifest} | tr "\t" "\n" | grep -n ${sampleid} | cut -d ":" -f1`
cut -f$COL ${manifest} | sed '1d' > tmp_sampleID

#get R1 file path from manifest file
COL=`head -n1 ${manifest} | tr "\t" "\n" | grep -n ${R1_files} | cut -d ":" -f1`
cut -f$COL ${manifest} | sed '1d' > tmp_R1

#get R2 file path from manifest file
[ ! ${singleEnd} ] && COL=`head -n1 ${manifest} | tr "\t" "\n" | grep -n ${R2_files} | cut -d ":" -f1`
[ ! ${singleEnd} ] && cut -f$COL ${manifest} | sed '1d' > tmp_R2

#get barcode
COL=`head -n1 ${metadata} | tr "\t" "\n" | grep -n ${barcode} | cut -d ":" -f1`
cut -f$COL ${metadata} | sed '1d' > tmp_barcodes

for i in $(cat tmp_R1) ; do echo $(zcat ${i}|wc -l)/4|bc >> tmp_reads ; done
for i in $(cat tmp_R1) ; do zcat ${i} > reads; head -n 1 reads | cut -d ':' -f1 >> tmp_sequencer; rm reads;  done

# Count and summary
paste -d ' ' tmp_barcodes tmp_R1 | sed "s/^/zgrep -c ':N:0:/g" | sed "s/ /' /3" | sh > tmp_barcodes_count_R1
[ ! ${singleEnd} ] && paste -d ' ' tmp_barcodes tmp_R2 | sed "s/^/zgrep -c ':N:0:/g" | sed "s/ /' /3" | sh > tmp_barcodes_count_R2
paste -d ' ' tmp_sequencer tmp_R1 | sed "s/^/zgrep -c '/g" | sed "s/ /' /3" | sh > tmp_sequencer_count_R1
[ ! ${singleEnd} ] && paste -d ' ' tmp_sequencer tmp_R2 | sed "s/^/zgrep -c '/g" | sed "s/ /' /3" | sh > tmp_sequencer_count_R2
for i in $(cat tmp_R1) ; do zgrep -c "${primerF}" ${i} >>  tmp_primer_count_R1 ; done
[ ! ${singleEnd} ] && for i in $(cat tmp_R2) ; do zgrep -c "${primerR}" ${i} >>  tmp_primer_count_R2 ; done
[ ! ${singleEnd} ] && paste -d '      ' tmp_sampleID tmp_barcodes tmp_reads tmp_barcodes_count_R1 tmp_barcodes_count_R2 tmp_sequencer_count_R1 tmp_sequencer_count_R2 tmp_primer_count_R1 tmp_primer_count_R2 > tmp_output
[ ${singleEnd} ] && paste -d '	' tmp_sampleID tmp_barcodes tmp_reads tmp_barcodes_count_R1 tmp_sequencer_count_R1 tmp_primer_count_R1 > tmp_output

# Barcode verification
valcol=4 #Count_in_R1
[ ${singleEnd} ] && addcol=7 #column = Perc_correct_barecode_R1
[ ! ${singleEnd} ] && addcol=10 #column = Perc_correct_barecode_R1
awk -v OFS='\t' -v "newcol=${addcol}" -v "col=${valcol}" '{$newcol = ($col != 0) ? sprintf("%.0f", $col/$3*100) : "0"}1' tmp_output | sort -nr -k ${addcol} > tmp_output.tmp && mv tmp_output.tmp tmp_output
if [ $(tail -n 1 tmp_output | awk -F "\t" -v "newcol=${addcol}" '{print $newcol}') -ge ${barcode_filter} ] ; then touch barcodes_R1.ok ; else touch barcodes_R1.bad ; fi
[ ! ${singleEnd} ] && awk -v OFS='\t' -v "newcol=$(($addcol+1))" -v "col=$(($valcol+1))" '{$newcol = ($col != 0) ? sprintf("%.0f", $col/$3*100) : "0"}1' tmp_output | sort -nr -k $(($addcol+1)) > tmp_output.tmp && mv tmp_output.tmp tmp_output
[ ! ${singleEnd} ] && if [ $(tail -n 1 tmp_output | awk -F"\t" -v "newcol=$(($addcol+1))" '{print $newcol}') -ge ${barcode_filter} ] ; then touch barcodes_R2.ok ; else touch barcodes_R2.bad ; fi

# Unique sequencer verification (Sequencer_R1, Sequencer_R2)
[ ${singleEnd} ] && valcol=5;addcol=8 #column = Perc_uniq_sequencer_R1 
[ ! ${singleEnd} ] && valcol=6;addcol=12 #column = Perc_uniq_sequencer_R1
awk -v OFS='\t' -v "newcol=${addcol}" -v "col=${valcol}" '{$newcol = ($col != 0) ? sprintf("%.0f", $col/$3*100) : "0"}1' tmp_output | sort -nr -k ${addcol} > tmp_output.tmp && mv tmp_output.tmp tmp_output
if [ $(tail -n 1 tmp_output | awk -F "\t" -v "newcol=${addcol}" '{print $newcol}') -ge 100 ] ; then touch sequencer_R1.ok ; else touch sequencer_R1.bad ; fi
[ ! ${singleEnd} ] && awk -v OFS='\t' -v "newcol=$(($addcol+1))" -v "col=$(($valcol+1))" '{$newcol = ($col != 0) ? sprintf("%.0f", $col/$3*100) : "0"}1' tmp_output | sort -nr -k $(($addcol+1)) > tmp_output.tmp && mv tmp_output.tmp tmp_output
[ ! ${singleEnd} ] && if [ $(tail -n 1 tmp_output | awk -F "\t" -v "newcol=$(($addcol+1))" '{print $newcol}') -ge 100 ] ; then touch sequencer_R2.ok ; else touch sequencer_R2.bad ; fi

# Primer verification (PrimerF_in_R1)
[ ${singleEnd} ] && valcol=6;addcol=9 #column = Perc_primer_R1
[ ! ${singleEnd} ] && valcol=8;addcol=14 #column = Perc_primer_R1
awk -v OFS='\t' -v "newcol=${addcol}" -v "col=${valcol}" '{$newcol = ($col != 0) ? sprintf("%.0f", $col/$3*100) : "0"}1' tmp_output | sort -nr -k ${addcol} > tmp_output.tmp && mv tmp_output.tmp tmp_output
if [ $(tail -n 1 tmp_output | awk -F "\t" -v "newcol=${addcol}" '{print $newcol}') -ge ${primer_filter} ] ; then touch primer_R1.ok ; else touch primer_R1.bad ; fi
[ ! ${singleEnd} ] && awk -v OFS='\t' -v "newcol=$(($addcol+1))" -v "col=$(($valcol+1))" '{$newcol = ($col != 0) ? sprintf("%.0f", $col/$3*100) : "0"}1' tmp_output | sort -nr -k $(($addcol+1)) > tmp_output.tmp && mv tmp_output.tmp tmp_output
[ ! ${singleEnd} ] && if [ $(tail -n 1 tmp_output | awk -F "\t" -v "newcol=$(($addcol+1))" '{print $newcol}') -ge ${primer_filter} ] ; then touch primer_R2.ok ; else touch primer_R2.bad ; fi

# Results saving and remove all temporary files
[ ${singleEnd} ] && if [[ -e barcodes_R1.ok && sequencer_R1.ok && -e primer_R1.ok ]] ; then touch ${verif_ok} ; else touch ${verif_bad} && echo '@@@ --- All verifications are not satisfied --- @@@' ; fi
[ ! ${singleEnd} ] && if [[ -e barcodes_R1.ok && -e barcodes_R2.ok && -e sequencer_R1.ok && -e sequencer_R2.ok && -e primer_R1.ok && -e primer_R2.ok ]] ; then touch ${verif_ok} ; else touch ${verif_bad} && echo '@@@ --- All verifications are not satisfied --- @@@' ; fi
mv tmp_output ${summary}
[ ${singleEnd} ] && sed -i '1 i\SampleID\tBarcode\tReads_count\tCount_in_R1\tSequencer_R1\tPrimerF_in_R1\tPerc_correct_barcode_R1\tPerc_uniq_sequencer_R1\tPerc_primerF_R1' ${summary}
[ ! ${singleEnd} ] && sed -i '1 i\SampleID\tBarcode\tReads_count\tCount_in_R1\tCount_in_R2\tSequencer_R1\tSequencer_R2\tPrimerF_in_R1\tPrimerR_in_R2\tPerc_correct_barcode_R1\tPerc_correct_barcode_R2\tPerc_uniq_sequencer_R1\tPerc_uniq_sequencer_R2\tPerc_primerF_R1\tPerc_primerR_R2' ${summary}
rm tmp*
