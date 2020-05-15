#!/usr/bin/env bash
###############################################################################
##                                                                           ##
## Purpose of script: Check integrity of raw data                            ##
##                                                                           ##
###############################################################################

args=("$@")
manifest=${args[0]}
metadata=${args[1]}
#Replace degenerated bases by dots in primers
primerF=$(echo "${args[2]}"|sed -e 's/[RYSWKMBDHVN]/\./g')
primerR=$(echo "${args[3]}"|sed -e 's/[RYSWKMBDHVN]/\./g')
summary=${args[4]}
barcode=${args[5]}
sampleid=${args[6]}

singleEnd=${args[12]}
if ${singleEnd}; then
    R1_files=${args[7]}
else 
    R1_files=${args[8]}
    R2_files=${args[9]}
fi

barcode_filter=${args[10]}
primer_filter=${args[11]}

verif_bad="data_integrity_control.bad"

# Verify if metadata file contains NA values
grep -P "NA\t|\tNA" ${metadata}
[ $? -eq 0 ] && { 
   echo "NA values found in $metadata, please remove them before running SAMBA"; 
   touch ${verif_bad} && echo "@@@ --- All controls are not satisfied --- @@@"; 
   exit 0; 
}

# Creation of temporary files
#get sample id list from manifest file
COL=$(head -n1 ${manifest}|tr "\t" "\n"|grep -n ${sampleid}|cut -d ":" -f1)
cut -f$COL ${manifest}|sed '1d' > tmp_sampleID

#get R1 file path from manifest file
COL=$(head -n1 ${manifest} | tr "\t" "\n" | grep -n ${R1_files} | cut -d ":" -f1)
cut -f$COL ${manifest} | sed '1d' > tmp_R1

#get R2 file path from manifest file
if ! ${singleEnd}; then
   COL=`head -n1 ${manifest} | tr "\t" "\n" | grep -n ${R2_files} | cut -d ":" -f1`
   cut -f$COL ${manifest} | sed '1d' > tmp_R2
fi
#get barcode
COL=$(head -n1 ${metadata} | tr "\t" "\n" | grep -n ${barcode} | cut -d ":" -f1)
cut -f$COL ${metadata} | sed '1d' > tmp_barcodes

for i in $(cat tmp_R1) ; do echo $(zcat ${i}|wc -l)/4|bc >> tmp_reads ; done
for i in $(cat tmp_R1) ; do zcat ${i} > reads; head -n 1 reads | cut -d ':' -f1 >> tmp_sequencer; rm reads;  done

# Count and summary
paste -d ' ' tmp_barcodes tmp_R1 | sed "s/^/zgrep -c ':N:0:/g" | sed "s/ /' /3" | sh > tmp_barcodes_count_R1
if ! ${singleEnd}; then
   paste -d ' ' tmp_barcodes tmp_R2 | sed "s/^/zgrep -c ':N:0:/g" | sed "s/ /' /3" | sh > tmp_barcodes_count_R2
fi
paste -d ' ' tmp_sequencer tmp_R1 | sed "s/^/zgrep -c '/g" | sed "s/ /' /3" | sh > tmp_sequencer_count_R1
if ! ${singleEnd}; then
   paste -d ' ' tmp_sequencer tmp_R2 | sed "s/^/zgrep -c '/g" | sed "s/ /' /3" | sh > tmp_sequencer_count_R2
fi

for i in $(cat tmp_R1) ; do zgrep -c "${primerF}" ${i} >>  tmp_primer_count_R1 ; done
if ! ${singleEnd}; then
   for i in $(cat tmp_R2) ; do zgrep -c "${primerR}" ${i} >>  tmp_primer_count_R2 ; done
   paste -d '      ' tmp_sampleID tmp_barcodes tmp_reads tmp_barcodes_count_R1 tmp_barcodes_count_R2 tmp_sequencer_count_R1 tmp_sequencer_count_R2 tmp_primer_count_R1 tmp_primer_count_R2 > tmp_output
else
   paste -d '	' tmp_sampleID tmp_barcodes tmp_reads tmp_barcodes_count_R1 tmp_sequencer_count_R1 tmp_primer_count_R1 > tmp_output
fi

# Barcode verification
valcol=4 #Count_in_R1
if ${singleEnd}; then
   addcol=7 #column = Perc_correct_barecode_R1
else
    addcol=10 #column = Perc_correct_barecode_R1
fi
awk -v OFS='\t' -v "newcol=${addcol}" -v "col=${valcol}" '{$newcol = ($col != 0) ? sprintf("%.0f", $col/$3*100) : "0"}1' tmp_output | sort -nr -k ${addcol} > tmp_output.tmp && mv tmp_output.tmp tmp_output
if [ $(tail -n 1 tmp_output | awk -F "\t" -v "newcol=${addcol}" '{print $newcol}') -ge ${barcode_filter} ] ; then touch barcodes_R1.ok ; else touch barcodes_R1.bad ; fi
if ! ${singleEnd}; then
   awk -v OFS='\t' -v "newcol=$(($addcol+1))" -v "col=$(($valcol+1))" '{$newcol = ($col != 0) ? sprintf("%.0f", $col/$3*100) : "0"}1' tmp_output | sort -nr -k $(($addcol+1)) > tmp_output.tmp && mv tmp_output.tmp tmp_output
   if [ $(tail -n 1 tmp_output | awk -F"\t" -v "newcol=$(($addcol+1))" '{print $newcol}') -ge ${barcode_filter} ] ; then 
      touch barcodes_R2.ok 
   else 
      touch barcodes_R2.bad
   fi
fi

# Unique sequencer verification (Sequencer_R1, Sequencer_R2)
if ${singleEnd}; then 
   valcol=5;addcol=8 #column = Perc_uniq_sequencer_R1 
else
   valcol=6;addcol=12 #column = Perc_uniq_sequencer_R1
fi
awk -v OFS='\t' -v "newcol=${addcol}" -v "col=${valcol}" '{$newcol = ($col != 0) ? sprintf("%.0f", $col/$3*100) : "0"}1' tmp_output | sort -nr -k ${addcol} > tmp_output.tmp && mv tmp_output.tmp tmp_output
if [ $(tail -n 1 tmp_output | awk -F "\t" -v "newcol=${addcol}" '{print $newcol}') -ge 100 ] ; then touch sequencer_R1.ok ; else touch sequencer_R1.bad ; fi
if ! ${singleEnd}; then
   awk -v OFS='\t' -v "newcol=$(($addcol+1))" -v "col=$(($valcol+1))" '{$newcol = ($col != 0) ? sprintf("%.0f", $col/$3*100) : "0"}1' tmp_output | sort -nr -k $(($addcol+1)) > tmp_output.tmp && mv tmp_output.tmp tmp_output
   if [ $(tail -n 1 tmp_output | awk -F "\t" -v "newcol=$(($addcol+1))" '{print $newcol}') -ge 100 ] ; then
      touch sequencer_R2.ok
   else 
      touch sequencer_R2.bad
   fi
fi

# Primer verification (PrimerF_in_R1)
if ${singleEnd}; then
   valcol=6
   addcol=9 #column = Perc_primer_R1
else 
   valcol=8
   addcol=14 #column = Perc_primer_R1
fi
awk -v OFS='\t' -v "newcol=${addcol}" -v "col=${valcol}" '{$newcol = ($col != 0) ? sprintf("%.0f", $col/$3*100) : "0"}1' tmp_output | sort -nr -k ${addcol} > tmp_output.tmp && mv tmp_output.tmp tmp_output
if [ $(tail -n 1 tmp_output | awk -F "\t" -v "newcol=${addcol}" '{print $newcol}') -ge ${primer_filter} ] ; then touch primer_R1.ok ; else touch primer_R1.bad ; fi
if ! ${singleEnd}; then
   awk -v OFS='\t' -v "newcol=$(($addcol+1))" -v "col=$(($valcol+1))" '{$newcol = ($col != 0) ? sprintf("%.0f", $col/$3*100) : "0"}1' tmp_output | sort -nr -k $(($addcol+1)) > tmp_output.tmp && mv tmp_output.tmp tmp_output
   if [ $(tail -n 1 tmp_output | awk -F "\t" -v "newcol=$(($addcol+1))" '{print $newcol}') -ge ${primer_filter} ] ; then 
      touch primer_R2.ok
   else 
      touch primer_R2.bad
   fi
fi

# Results saving and remove all temporary files
if ${singleEnd}; then
   if [[ -e barcodes_R1.ok && -e sequencer_R1.ok && -e primer_R1.ok ]] ; then 
      echo '@@@ --- Your data passed the control --- @@@'
   else 
      touch ${verif_bad} && echo '@@@ --- All controls are not satisfied --- @@@'
   fi
else
   if [[ -e barcodes_R1.ok && -e barcodes_R2.ok && -e sequencer_R1.ok && -e sequencer_R2.ok && -e primer_R1.ok && -e primer_R2.ok ]] ; then 
      echo '@@@ --- Your data passed the control --- @@@' 
   else 
      touch ${verif_bad} && echo '@@@ --- All controls are not satisfied --- @@@'
   fi
fi
mv tmp_output ${summary}
if ${singleEnd}; then
   sed -i '1 i\SampleID\tBarcode\tReads_count\tCount_in_R1\tSequencer_R1\tPrimerF_in_R1\tPerc_correct_barcode_R1\tPerc_uniq_sequencer_R1\tPerc_primerF_R1' ${summary}
else
   sed -i '1 i\SampleID\tBarcode\tReads_count\tCount_in_R1\tCount_in_R2\tSequencer_R1\tSequencer_R2\tPrimerF_in_R1\tPrimerR_in_R2\tPerc_correct_barcode_R1\tPerc_correct_barcode_R2\tPerc_uniq_sequencer_R1\tPerc_uniq_sequencer_R2\tPerc_primerF_R1\tPerc_primerR_R2' ${summary}
fi

rm tmp*
