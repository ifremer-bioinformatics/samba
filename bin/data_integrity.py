#!/usr/bin/env python3

from __future__ import print_function
import re
import gzip
import argparse
import sys
from Bio import SeqIO
import multiprocessing as mp

# For errors / warnings
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def getArgs():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-a',dest="manifest",type=str,required=True,help='q2_manifest file')
    parser.add_argument('-e',dest="metadata",type=str,required=True,help='q2_metadata file')
    parser.add_argument('-p',dest="primer",type=int,default=70,help='Percentage of primers supposed to be found in raw reads [%(default)s]')
    parser.add_argument('-c',dest="control",type=str,help='Delimited list of control (comma separated)')
    parser.add_argument('-r',dest="readtype",type=str,default='paired', choices=['paired','single'],help='Sequencing format [%(default)s]')
    parser.add_argument('-t',dest="threads",type=int,default=4,help='Number of threads [%(default)s]')

    arg = parser.parse_args()

    return arg

def collect_check_metadata(metadata_file, data_type):
    collect_data = {}
    # read metadata line by line
    for l in open(metadata_file, 'r'):
        # ignore header
        if not l.startswith('sampleid'):
            sampleid = l.split()[0]
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
                sampleid, R1, R2 = l.split()
                try:
                    collect_data[sampleid]['R1'] = R1
                    collect_data[sampleid]['R2'] = R2
                except KeyError:
                    eprint('ERROR: '+sampleid+' from q2_manifest is absent in q2_metadata')
                    exit(1)
            else:
                sampleid, R1 = l.split()
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
    collect_data, data_type = collect_check_metadata(args.metadata, args.readtype)
    # Collect reads location from manifest
    print("Step 2 - parse and collect data from q2_manifest")
    collect_data = collect_manifest(args.manifest, collect_data, data_type)
    # Check fastq integrity
    print("Step 3 - check fastq integrity")
    pool = mp.Pool(args.threads)
    collect_data_para = pool.starmap(check_fastq, [(collect_data, sample,data_type) for sample in collect_data.keys()])
    pool.close()
    # Clean all this mess
    collect_data_out = {}
    for result in collect_data_para:
        for sample, vals in result.items():
            collect_data_out[sample] = vals
    # Report fastq intergrity before validation
    # Allow exploration by user for a better understanding
    print("Step 4 - write integrity report")
    write_report(collect_data_out, data_type)
    # Validate manifest, metadata and fastq
    print("Step 5 - validation of data")
    controls = [str(control) for control in args.control.split(',')]
    print(controls)
    for sample, val in collect_data_out.items():
        # print(sample, val)
        # check samples
        if val['reads_count_R1'] != val['reads_count_R2']:
            eprint('ERROR: '+sampleid+' from q2_manifest is absent in q2_metadata')
            exit(1)

            
        if not sample in controls:

            # if val['nb_instrument_R1'] != val['reads_count_R2']:
            if val['perc_primer_R1'] <= args.primer or val['perc_primer_R2'] <= args.primer:
                eprint('ERROR: '+sampleid+' from q2_manifest is absent in q2_metadata')
                exit(1)


    # for k, v in collect_data_out.items()
        # -> primer filter
        # -> primer filter expcet for controls
        # -> barcode filter: should be == TRUE
        # -> sequencer name filter: should == 1

    # print("Step 6 - rewrite manifest and metadata (to sort and remove space)")

if __name__ == '__main__':
    args = getArgs()
    main(args)
