#!/usr/bin/env python3

from __future__ import print_function
import re
import os
import sys
import gzip
import argparse
import subprocess
import multiprocessing as mp
from Bio import SeqIO

# For errors / warnings
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def getArgs():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-a',dest="manifest",type=str,required=True,help='q2_manifest file')
    parser.add_argument('-e',dest="metadata",type=str,required=True,help='q2_metadata file')
    parser.add_argument('-p',dest="primer",type=int,default=70,help='Percentage of primers supposed to be found in raw reads [%(default)s]')
    parser.add_argument('-c',dest="control",type=str, default='',help='Delimited list of control (comma separated)')
    parser.add_argument('-r',dest="readtype",type=str,default='paired',choices=['paired','single'],help='Sequencing format [%(default)s]')
    parser.add_argument('-t',dest="threads",type=int,default=4,help='Number of threads [%(default)s]')

    arg = parser.parse_args()

    return arg

def collect_check_metadata(metadata_file, data_type):
    collect_data = {}
    # read metadata line by line
    for l in open(metadata_file, 'r'):
        # ignore header
        if not l.startswith('sampleid'):
            sampleid = re.split(r'\t', l.rstrip('\n'))[0]
            collect_data[sampleid] = {}
            collect_data[sampleid]['barcode'] = re.split(r'\t', l.rstrip('\n'))[1]
            if data_type == 'paired':
                collect_data[sampleid]['primerF'] = re.split(r'\t', l.rstrip('\n'))[2]
                collect_data[sampleid]['primerR'] = re.split(r'\t', l.rstrip('\n'))[3]
                collect_data[sampleid]['vars'] = re.split(r'\t', l.rstrip('\n'))[4:]
            else:
                collect_data[sampleid]['primerF'] = re.split(r'\t', l.rstrip('\n'))[2]
                collect_data[sampleid]['vars'] = re.split(r'\t', l.rstrip('\n'))[3:]
            # check that var(s) didn't contains any NA
            for var in collect_data[sampleid]['vars']:
                if var == 'NA' or var == '':
                    eprint('ERROR: '+sampleid+' has NA value(s) in q2_metadata, please remove them before running SAMBA')
                    exit(1)
    # return results
    return collect_data, data_type

def collect_manifest(manifest_file, collect_data, data_type):
    # in order to check that metadata and manifest have the same number of lines
    metadata_size = len(collect_data)
    manifest_size = 0
    # read manifest
    for l in open(manifest_file, 'r'):
        # ignore header
        if not l.startswith('sample-id'):
            if data_type == 'paired':
                try:
                    sampleid, R1, R2 = re.split(r'\t', l.rstrip('\n'))
                except ValueError:
                    size = len(re.split(r'\t', l.rstrip('\n')))
                    eprint('ERROR: q2_manifest contains '+str(size)+ 'column(s) instead of 3')
                    exit(1)

                # check that file path are good
                check_fastq_path(R1)
                check_fastq_path(R2)
                # add to collect_data
                try:
                    collect_data[sampleid]['R1'] = R1
                    collect_data[sampleid]['R2'] = R2
                except KeyError:
                    eprint('ERROR: '+sampleid+' from q2_manifest is absent in q2_metadata')
                    exit(1)
            else:
                try:
                    sampleid, R1 = l.split()
                except ValueError:
                    size = len(re.split(r'\t', l.rstrip('\n')))
                    eprint('ERROR: q2_manifest contains '+str(size)+ 'column(s) instead of 2')
                    exit(1)
                # check that file path are good
                check_fastq_path(R1)
                # add to collect_data
                try:
                    collect_data[sampleid]['R1'] = R1
                except KeyError:
                    eprint('ERROR: '+sampleid+' from q2_manifest is absent in q2_metadata')
                    exit(1)
            # increase line counter
            manifest_size += 1
    # check that manifest have the same size as metadata
    if manifest_size != metadata_size:
        eprint('ERROR: q2_manifest and q2_metadata did not have the same number of lines')
        exit(1)
    # return update dict
    return collect_data

def check_fastq(collect_data, sample, data_type):
    out = {}
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
        # chech that barcode are the same as in metadata
        if len(index) == 1:
            if index[0] == barcode:
                collect_data[sample]['barcode_seq_'+r]  = 'TRUE'
    # return updated dict
    out[sample] = collect_data[sample]
    return out

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

def check_fastq_path(path):
    if not os.path.isfile(path):
        eprint('ERROR: ' + path + ' from q2_manifest not exit. Wrong path?')
        exit(1)

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
            report.write(sample+'\t'+'{reads_count_R1}\t{reads_count_R2}\t{barcode}\t{nb_barcode_R1}\t{barcode_seq_R1}\t{nb_barcode_R2}\t{barcode_seq_R2}\t{nb_instrument_R1}\t{nb_instrument_R2}\t{primer_R1}\t{perc_primer_R1}\t{primer_R2}\t{perc_primer_R2}\n'.format(**val))
        else:
            report.write(sample+'\t'+'{reads_count_R1}\t{barcode}\t{nb_barcode_R1}\t{barcode_seq_R1}\t{nb_instrument_R1}\t{primer_R1}\t{perc_primer_R1}\n'.format(**val))

def integrity_validation(collect_data_out,control_list, data_type, primer_threshold):
    # collect control(s)
    controls = [str(control) for control in control_list.split(',')]
    # let's check the integrity
    for sample_id, val in collect_data_out.items():
        if data_type == 'paired':
            # 1 - check reads count R1 vs R2
            if val['reads_count_R1'] != val['reads_count_R2']:
                eprint('ERROR: different number of reads between R1 and R2 for ' + sample_id)
                exit(1)
            # 2 - check single sequencing instrument in reads
            if val['nb_instrument_R1'] != 1 or val['nb_instrument_R2'] != 1:
                eprint('ERROR: multiple sequencing machines detected in ' + sample_id)
                exit(1)
            # 3 - check primer percentage found, except for control(s)
            if not sample_id in controls:
                if val['perc_primer_R1'] < primer_threshold or val['perc_primer_R2'] < primer_threshold:
                    eprint('ERROR: ' + sample_id + " did not reach the minimum threshold for primer percentage [" + str(primer_threshold) +"%]")
                    exit(1)
        else:
            # 1 - check single sequencing instrument in reads
            if val['nb_instrument_R1'] != 1:
                eprint('ERROR: multiple sequencing machines detected in ' + sample_id)
                exit(1)
            # 2 - check primer percentage found, except for control(s)
            if not sample_id in controls:
                if val['perc_primer_R1'] < primer_threshold:
                    eprint('ERROR: ' + sample_id + " did not reach the minimum threshold for primer percentage [" + str(primer_threshold) + "%]")
                    exit(1)

def sort_process(input_file, output_file):
    cmd_sort = '(head -n 1 {in_file} && tail -n +2 {in_file} | sort) > {out_file}'.format(in_file=input_file, out_file=output_file)
    p = subprocess.Popen(cmd_sort, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = p.communicate()
    if p.returncode != 0:
        eprint('ERROR: sort metadata/manifest failed')
        raise Exception(stderr)
        exit(1)

def main(args):

    # 1 - Collect metadata infos and check variable names
    print("Step 1 - parse and collect data from q2_metadata")
    collect_data, data_type = collect_check_metadata(args.metadata, args.readtype)

    # 2 - Collect reads location from manifest
    print("Step 2 - parse and collect data from q2_manifest")
    collect_data = collect_manifest(args.manifest, collect_data, data_type)

    # 3 - Check fastq integrity
    print("Step 3 - collect fastq integrity")
    pool = mp.Pool(args.threads)
    collect_data_para = pool.starmap(check_fastq, [(collect_data, sample,data_type) for sample in collect_data.keys()])
    pool.close()
    # Clean all this mess
    collect_data_out = {}
    for result in collect_data_para:
        for sample, vals in result.items():
            collect_data_out[sample] = vals

    # 4 - Report fastq intergrity before validation
    # Allow exploration by user for a better understanding
    print("Step 4 - write integrity report")
    write_report(collect_data_out, data_type)

    # 5 - Validate manifest, metadata and fastq
    print("Step 5 - validation of data")
    integrity_validation(collect_data_out, args.control, data_type, args.primer)

    # 6 - Sort manifest and metadata
    print("Step 6 - Sort manifest and metadata by sample id")
    print("\tsort manifest...")
    output_file_manifest = args.manifest + ".sort"
    sort_process(args.manifest, output_file_manifest)
    print("\tsort metadata...")
    output_file_metadata = args.metadata + ".sort"
    sort_process(args.metadata, output_file_metadata)

if __name__ == '__main__':
    args = getArgs()
    main(args)
