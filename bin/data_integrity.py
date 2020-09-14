#!/usr/bin/env python3

from __future__ import print_function
# import os
# import glob
import re
import gzip
import argparse
import sys
from Bio import SeqIO
from joblib import Parallel, delayed

# For errors / warnings
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def getArgs():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-a',dest="manifest",type=str,required=True,help='q2_manifest file')
    parser.add_argument('-e',dest="metadata",type=str,required=True,help='q2_metadata file')
    # parser.add_argument('-s',dest="summary",type=str,required=True,help='Summary file')
    # parser.add_argument('-b',dest="barcode",type=str,required=True,help='Barcode sequence')
    # parser.add_argument('-i',dest="sampleid",type=str,required=True,help='Sample ID')
    parser.add_argument('-t',dest="type",type=str,default='paired', choices=['paired','single'],help='Sequencing format [%(default)s]')
    # parser.add_argument('-x',dest="barcodefilter",type=int,required=True,help='Barcode percentage filter')
    # parser.add_argument('-y',dest="primerfilter",type=int,required=True,help='Primer percentage filter')

    arg = parser.parse_args()
    return arg

def collect_check_metadata(metadata_file, data_type):
    collect_data = {}
    # read metadata line by line
    for l in open(metadata_file, 'r'):
        sampleid = l.split()[0]
        # ignore header
        if sampleid != 'sampleid':
            collect_data[sampleid] = {}
            collect_data[sampleid]['barcode'] = l.split()[1]
            if data_type == 'paired':
                collect_data[sampleid]['primerF'] = l.split()[2]
                collect_data[sampleid]['primerR'] = l.split()[3]
                collect_data[sampleid]['vars'] = l.split()[4:]
            else:
                collect_data[sampleid]['primerF'] = l.split()[2]
                collect_data[sampleid]['vars'] = l.split()[3:]
            # check that var(s) didn't contains any NA
            for var in collect_data[sampleid]['vars']:
                if var == 'NA':
                    eprint('ERROR: '+sampleid+' has NA value(s) in $metadata, please remove them before running SAMBA')
                    verify_bad=open('verify.bad', 'w')
                    verify_bad.close()
                    exit(1)
    # process file validator
    verify_ok=open('verify.ok', 'w')
    verify_ok.close()
    # return results
    return collect_data, data_type

def collect_manifest(manifest_file, collect_data, data_type):
    for l in open(manifest_file, 'r'):
        # ignore header
        if not l.startswith('sample-id'):
            if data_type == 'paired':
                sampleid, R1, R2 = l.split()
                collect_data[sampleid]['R1'] = R1
                collect_data[sampleid]['R2'] = R2
            else:
                sampleid, R1 = l.split()
                collect_data[sampleid]['R1'] = R1
    return collect_data

def check_fastq(collect_data, data_type):
    for sample in collect_data.keys():
        print("\tanalyse sample: "+sample)
        barcode = collect_data[sample]['barcode']
        # by default, single-end
        reads = ['R1']
        # in case of paired, add R2 analysis
        if data_type == 'paired':
            reads.append('R2')
        # check fastq
        for r in reads:
            R=collect_data[sample][r]
            if r == 'R1':
                primer=collect_data[sample]['primerF']
            else:
                primer=collect_data[sample]['primerR']
            instrument, index, reads_count, primer = read_fastq(R, primer)
            collect_data[sample]['reads_count_'+r] = reads_count
            collect_data[sample]['instrument_'+r] = instrument
            collect_data[sample]['nb_instrument_'+r] = len(instrument)
            collect_data[sample]['index_'+r] = index
            collect_data[sample]['primer_'+r] = primer
            collect_data[sample]['perc_primer_'+r] = round(primer * 100 / reads_count, 2)
            collect_data[sample]['nb_barcode_'+r] = len(index)
            collect_data[sample]['barcode_seq_'+r] = 'FALSE'
            if len(index) == 1:
                if index[0] == barcode:
                    collect_data[sample]['barcode_seq_'+r]  = 'TRUE'
    # return updated dict
    return collect_data

def read_fastq(fastq, primer):
    primer = re.sub(r"([RYSWKMBDHVN])", r".", primer)
    reads_count = 0
    primers_count = 0
    instrument = []
    index = []
    with gzip.open(fastq, "rt") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                # check instrument name
                h_1, h_2 = record.description.split()
                instrument_name = h_1.split(':')[0]
                if not instrument_name in instrument:
                    instrument.append(instrument_name)
                # check index
                sequence_index = h_2.split(':')[3]
                if not sequence_index in index:
                    index.append(sequence_index)
                # check primers
                if re.search(primer, str(record.seq)):
                    primers_count += 1
                # increase read count
                reads_count += 1
    return instrument, index, reads_count, primers_count

def write_report(collect_data, data_type):
    # open output file for writting
    report = open('data_integrity.txt', 'w')
    # header
    header_S = ['SampleID', 'Reads_R1', 'Barcode', 'Uniq_in_R1', 'Same_as_ref_R1', 'Uniq_sequencer_R1', 'PrimerF_in_R1', 'Perc_primerF_R1']
    header_P = ['SampleID', 'Reads_R1', 'Reads_R2', 'Barcode', 'Uniq_in_R1', 'Same_as_ref_R1', 'Uniq_in_R2', 'Same_as_ref_R2', 'Uniq_sequencer_R1', 'Uniq_sequencer_R2', 'PrimerF_in_R1', 'Perc_primerF_R1', 'PrimerR_in_R2', 'Perc_primerR_R2']
    if data_type == 'paired':
        report.write('\t'.join(header_P)+'\n')
    else:
        report.write('\t'.join(header_S)+'\n')
    for sample, val in collect_data.items():
        if data_type == 'paired':
            header_P = ['SampleID', 'Reads_R1', 'Reads_R2', 'Barcode', 'Uniq_in_R1', 'Same_as_ref_R1', 'Uniq_in_R2', 'Same_as_ref_R2', 'Uniq_sequencer_R1', 'Uniq_sequencer_R2', 'PrimerF_in_R1', 'Perc_primerF_R1', 'PrimerR_in_R2', 'Perc_primerR_R2']
            report.write(sample+'\t'+'{reads_count_R1}\t{reads_count_R2}\t{barcode}\t{nb_barcode_R1}\t{barcode_seq_R1}\t{nb_barcode_R2}\t{barcode_seq_R2}\t{nb_instrument_R1}\t{nb_instrument_R2}\t{primer_R1}\t{perc_primer_R1}\t{primer_R2}\t{perc_primer_R2}\n'.format(**val))
        else:
            report.write(sample+'\t'+'{reads_count_R1}\t{barcode}\t{nb_barcode_R1}\t{barcode_seq_R1}\t{nb_instrument_R1}\t{primer_R1}\t{perc_primer_R1}\n'.format(**val))

def main(args):

    # Collect metadata infos and check variable names
    print("Step 1 - parse and collect data from q2_metadata")
    collect_data, data_type = collect_check_metadata(args.metadata, args.type)
    # Collect reads location from manifest
    print("Step 2 - parse and collect data from q2_manifest")
    collect_data = collect_manifest(args.manifest, collect_data, data_type)
    # Check fastq integrity
    print("Step 3 - check fastq integrity")
    collect_data = check_fastq(collect_data, data_type)
    # Validate and report fastq intergrity
    print("Step 4 - write integrity report")
    write_report(collect_data, data_type)

if __name__ == '__main__':
    args = getArgs()
    main(args)
