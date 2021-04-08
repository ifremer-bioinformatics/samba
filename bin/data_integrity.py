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
import xlrd

# For errors / warnings
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def getArgs():
    parser = argparse.ArgumentParser(description="")
    # parser.add_argument('-x', dest="xlsx", type=str, required=False, help='Data description (Excel file)')
    parser.add_argument('-a',dest="manifest",type=str,required=True,help='q2_manifest file')
    parser.add_argument('-e',dest="metadata",type=str,required=True,help='q2_metadata file')
    parser.add_argument('-p',dest="primer",type=int,default=70,help='Percentage of primers supposed to be found in raw reads [%(default)s]')
    parser.add_argument('-c',dest="control",type=str, default='',help='Delimited list of control (comma separated)')
    parser.add_argument('-r',dest="readtype",type=str,default='paired',choices=['paired','single'],help='Sequencing format [%(default)s]')
    parser.add_argument('-t',dest="threads",type=int,default=4,help='Number of threads [%(default)s]')

    arg = parser.parse_args()

    return arg

def import_metadata(metadata_file, data_type):
    samba_samples = {}
    # read metadata line by line
    for l in open(metadata_file, 'r'):
        # ignore header
        if not l.startswith('sampleid'):
            sampleid = re.split(r'\t', l.rstrip('\n'))[0]
            samba_samples[sampleid] = {}
            samba_samples[sampleid]['barcode'] = re.split(r'\t', l.rstrip('\n'))[1]
            samba_samples[sampleid]['primerF'] = re.split(r'\t', l.rstrip('\n'))[2]
            if data_type == 'paired':
                samba_samples[sampleid]['primerR'] = re.split(r'\t', l.rstrip('\n'))[3]
                samba_samples[sampleid]['vars'] = re.split(r'\t', l.rstrip('\n'))[4:]
            else:
                samba_samples[sampleid]['vars'] = re.split(r'\t', l.rstrip('\n'))[3:]
            # check that var(s) didn't contains any NA
            for var in samba_samples[sampleid]['vars']:
                if var == 'NA' or var == '':
                    eprint('ERROR: '+sampleid+' has NA value(s) in q2_metadata, please remove them before running SAMBA')
                    exit(1)
    # Overlap param for Cutadapt
    primers_lgth = [samba_samples[sampleid]['primerF']]
    if data_type == 'paired':
        primers_lgth = [samba_samples[sampleid]['primerF'], samba_samples[sampleid]['primerR']]
    cutadapt_overlap = len(min(primers_lgth, key=len)) - 1
    # return results
    return samba_samples, cutadapt_overlap

def import_manifest(manifest_file, metadata, data_type):
    # in order to check that metadata and manifest have the same number of lines
    metadata_size = len(metadata)
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
                    metadata[sampleid]['R1'] = R1
                    metadata[sampleid]['R2'] = R2
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
                    metadata[sampleid]['R1'] = R1
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
    return metadata

def check_fastq(samba_samples, sample, data_type):
    out = {}
    print("\tanalyse sample: "+sample)
    barcode = samba_samples[sample]['barcode']
    # by default, single-end
    reads = ['R1']
    # in case of paired, add R2 analysis
    if data_type == 'paired':
        reads.append('R2')
    # check fastq
    for r in reads:
        R=samba_samples[sample][r]
        if r == 'R1':
            primer=samba_samples[sample]['primerF']
        else:
            primer=samba_samples[sample]['primerR']
        instrument, index, reads_count, primer = read_fastq(R, primer)
        samba_samples[sample]['reads_count_' + r] = reads_count
        samba_samples[sample]['instrument_' + r] = instrument
        samba_samples[sample]['nb_instrument_' + r] = len(instrument)
        samba_samples[sample]['index_' + r] = index
        samba_samples[sample]['primer_' + r] = primer
        samba_samples[sample]['perc_primer_' + r] = round(primer * 100 / reads_count, 2)
        samba_samples[sample]['nb_barcode_' + r] = len(index)
        samba_samples[sample]['barcode_seq_' + r] = 'FALSE'
        # chech that barcode are the same as in metadata
        if len(index) == 1:
            if index[0] == barcode:
                samba_samples[sample]['barcode_seq_' + r]  = 'TRUE'
    # return updated dict
    out[sample] = samba_samples[sample]
    return out

def read_fastq(fastq, primer):
    primer = re.sub(r"([RYSWKMBDHVNI])", r".", primer)
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

def write_report(integrity_data, data_type):
    # open output file for writting
    report = open('data_integrity.txt', 'w')
    # header
    header_S = ['SampleID', 'Reads_R1', 'Barcode', 'Uniq_in_R1', 'Same_as_ref_R1', 'Uniq_sequencer_R1', 'PrimerF_in_R1', 'Perc_primerF_R1']
    header_P = ['SampleID', 'Reads_R1', 'Reads_R2', 'Barcode', 'Uniq_in_R1', 'Same_as_ref_R1', 'Uniq_in_R2', 'Same_as_ref_R2', 'Uniq_sequencer_R1', 'Uniq_sequencer_R2', 'PrimerF_in_R1', 'Perc_primerF_R1', 'PrimerR_in_R2', 'Perc_primerR_R2']
    if data_type == 'paired':
        report.write('\t'.join(header_P)+'\n')
    else:
        report.write('\t'.join(header_S)+'\n')
    for sample, val in integrity_data.items():
        if data_type == 'paired':
            report.write(sample+'\t'+'{reads_count_R1}\t{reads_count_R2}\t{barcode}\t{nb_barcode_R1}\t{barcode_seq_R1}\t{nb_barcode_R2}\t{barcode_seq_R2}\t{nb_instrument_R1}\t{nb_instrument_R2}\t{primer_R1}\t{perc_primer_R1}\t{primer_R2}\t{perc_primer_R2}\n'.format(**val))
        else:
            report.write(sample+'\t'+'{reads_count_R1}\t{barcode}\t{nb_barcode_R1}\t{barcode_seq_R1}\t{nb_instrument_R1}\t{primer_R1}\t{perc_primer_R1}\n'.format(**val))

def integrity_validation(integrity_data,control_list, data_type, primer_threshold):
    # collect control(s)
    controls = [str(control) for control in control_list.split(',')]
    # let's check the integrity
    for sample_id, val in integrity_data.items():
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
    metadata, cutadapt_overlap = import_metadata(args.metadata, args.readtype)

    # 2 - Collect reads location from manifest
    print("Step 2 - parse and collect data from q2_manifest")
    metadata = import_manifest(args.manifest, metadata, args.readtype)

    # 3 - Check fastq integrity
    print("Step 3 - check fastq integrity")
    pool = mp.Pool(args.threads)
    fastq = pool.starmap(check_fastq, [(metadata, sample, args.readtype) for sample in metadata.keys()])
    pool.close()
    # Clean all this mess
    integrity_data = {}
    for result in fastq:
        for sample, vals in result.items():
            integrity_data[sample] = vals

    # 4 - Report fastq intergrity before validation
    # Allow exploration by user for a better understanding
    print("Step 4 - write integrity report")
    write_report(integrity_data, args.readtype)

    # 5 - Validate manifest, metadata and fastq
    print("Step 5 - validation of data")
    integrity_validation(integrity_data, args.control, args.readtype, args.primer)

    # 6 - Sort manifest and metadata
    print("Step 6 - Sort manifest and metadata by sample id")
    print("\tsort manifest...")
    sorted_manifest = args.manifest + ".sort"
    sort_process(args.manifest, sorted_manifest)
    print("\tsort metadata...")
    sorted_metadata = args.metadata + ".sort"
    sort_process(args.metadata, sorted_metadata)

    # 7 - Output cutadapt overlap length
    cutovl = open('cutovl.txt','w')
    cutovl.write(str(cutadapt_overlap))
    cutovl.close()

if __name__ == '__main__':
    args = getArgs()
    main(args)
