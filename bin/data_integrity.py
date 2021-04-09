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
    parser.add_argument('-x', dest="xlsx", type=str, required=True, help='Data description (Excel file)')
    # parser.add_argument('-a',dest="manifest",type=str,required=True,help='q2_manifest file')
    # parser.add_argument('-e',dest="metadata",type=str,required=True,help='q2_metadata file')
    parser.add_argument('-p',dest="primer",type=int,default=70,help='Percentage of primers supposed to be found in raw reads [%(default)s]')
    parser.add_argument('-c',dest="control",type=str, default='',help='Delimited list of control (comma separated)')
    parser.add_argument('-r',dest="readtype",type=str,default='paired',choices=['paired','single'],help='Sequencing format [%(default)s]')
    parser.add_argument('-t',dest="threads",type=int,default=4,help='Number of threads [%(default)s]')

    arg = parser.parse_args()

    return arg

def import_metadata(metadata_file, data_type):
    samba_metadata = {}
    # read metadata line by line
    for l in open(metadata_file, 'r'):
        # ignore header
        if not l.startswith('sampleid'):
            sampleid = re.split(r'\t', l.rstrip('\n'))[0]
            samba_metadata[sampleid] = {}
            samba_metadata[sampleid]['barcode'] = re.split(r'\t', l.rstrip('\n'))[1]
            samba_metadata[sampleid]['primerF'] = re.split(r'\t', l.rstrip('\n'))[2]
            if data_type == 'paired':
                samba_metadata[sampleid]['primerR'] = re.split(r'\t', l.rstrip('\n'))[3]
                samba_metadata[sampleid]['vars'] = re.split(r'\t', l.rstrip('\n'))[4:]
            else:
                samba_metadata[sampleid]['vars'] = re.split(r'\t', l.rstrip('\n'))[3:]
            # check that var(s) didn't contains any NA
            for var in samba_metadata[sampleid]['vars']:
                if var == 'NA' or var == '':
                    eprint('ERROR: '+sampleid+' has NA value(s) in q2_metadata, please remove them before running SAMBA')
                    exit(1)
    # Overlap param for Cutadapt
    primers_lgth = [samba_metadata[sampleid]['primerF']]
    if data_type == 'paired':
        primers_lgth = [samba_metadata[sampleid]['primerF'], samba_metadata[sampleid]['primerR']]
    cutadapt_overlap = len(min(primers_lgth, key=len)) - 1
    # return results
    return samba_metadata, cutadapt_overlap

def import_manifest(manifest_file, data_type):
    samba_manifest = {}
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
                samba_manifest[sampleid] = {}
                samba_manifest[sampleid]['R1'] = R1
                samba_manifest[sampleid]['R2'] = R2
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
                samba_manifest[sampleid] = {}
                samba_manifest[sampleid]['R1'] = R1
    # return results
    return samba_manifest

def check_manifest_vs_metadata(manifest, metadata):
    # 1 - check samples names and number consistency
    manifest_samples = set()
    metadata_samples = set()
    # collect samples in manifest and metadata
    for k in manifest.keys():
        manifest_samples.add(k)
    for k in metadata.keys():
        metadata_samples.add(k)
    # stop if any difference
    if manifest_samples != metadata_samples:
        uniq_manifest = manifest_samples.difference(metadata_samples)
        uniq_metadata = metadata_samples.difference(manifest_samples)
        eprint('ERROR: Manifest and metadata did not have the same sample(s)')
        if len(uniq_manifest) > 0:
            eprint('ERROR: Sample(s) only found in manifest -> ', ','.join(uniq_manifest))
        if len(uniq_metadata) > 0:
            eprint('ERROR: Sample(s) only found in metadata -> ', ','.join(uniq_metadata))
        exit(1)
    # merge dict
    samba_samples = manifest.copy()
    for k in metadata.keys():
        for sk in metadata[k]:
            samba_samples[k][sk] = metadata[k][sk]
    # return merged
    return samba_samples

def check_fastq(samba_samples, sample, data_type):
    out = {}
    print("\t\tanalyse sample: "+sample)
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
        eprint('ERROR: ' + path + ' from q2_manifest do not exist. Wrong path?')
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

def import_excel(xlsx, single=False):
    samba_samples = {}
    # 1 - load data
    data = xlrd.open_workbook(xlsx)
    manifest = data.sheet_by_name('manifest')
    metadata = data.sheet_by_name('metadata')

    row_count_manifest = manifest.nrows
    col_count_manifest = manifest.ncols
    row_count_metadata = metadata.nrows

    # 2 - check samples names and number consistency
    manifest_samples = set()
    metadata_samples = set()
    # collect samples in manifest and metadata
    for row in range(1, row_count_manifest):
        manifest_samples.add(manifest.row_values(row, start_colx=0, end_colx=None)[0])
    for row in range(1, row_count_metadata):
        metadata_samples.add(metadata.row_values(row, start_colx=0, end_colx=None)[0])
    # stop if any difference
    if manifest_samples != metadata_samples:
        uniq_manifest = manifest_samples.difference(metadata_samples)
        uniq_metadata = metadata_samples.difference(manifest_samples)
        eprint('ERROR: Manifest and metadata did not have the same sample(s)')
        if len(uniq_manifest) > 0:
            eprint('ERROR: Sample(s) only found in manifest -> ', ','.join(uniq_manifest))
        if len(uniq_metadata) > 0:
            eprint('ERROR: Sample(s) only found in metadata -> ', ','.join(uniq_metadata))
        exit(1)

    # 3 - manifest
    ## 3.1 - check column number
    col = 3
    if single:
        col = 2
    if col_count_manifest != col:
        eprint('ERROR: manifest did not have the required number of columns')
        eprint('ERROR: ' + str(col_count_manifest) + ' instead of ' + str(col) + ' columns')
        exit(1)

    ## 3.2 - parse the manifest sheet
    for row in range(1, row_count_manifest):
        # read excel row by row
        if single:
            sample, R1 = manifest.row_values(row, start_colx=0, end_colx=None)
            # check that file path are good
            check_fastq_path(R1)
            # add to samba_samples
            samba_samples[sample] = {}
            samba_samples[sample]['R1'] = R1
        else:
            sample, R1, R2 = manifest.row_values(row, start_colx=0, end_colx=None)
            # check that file path are good
            check_fastq_path(R1)
            check_fastq_path(R2)
            # add to samba_samples
            samba_samples[sample] = {}
            samba_samples[sample]['R1'] = R1
            samba_samples[sample]['R2'] = R2

    # 4 - metadata
    ## 4.1 - parse the metadata sheet
    for row in range(1, row_count_metadata):
        tmp = metadata.row_values(row, start_colx=0, end_colx=None)
        sample = tmp[0]
        # store metadata
        samba_samples[sample]['barcode'] = tmp[1]
        samba_samples[sample]['primerF'] = tmp[2]
        if single:
            samba_samples[sample]['vars'] = tmp[3:]
        else:
            samba_samples[sample]['primerR'] = tmp[3]
            samba_samples[sample]['vars'] = tmp[4:]
        # check empty variable and replace \s
        for var in samba_samples[sample]['vars']:
            if var in ['NA', '', 'Nan', 'NaN', 'nan']:
                eprint(
                    'ERROR: ' + sample + ' has unset values (NA, NaN, empty cell...) in metadata, please fill them before running SAMBA')
                exit(1)

    ## 4.2 - get overlap size
    primers_lgth = [samba_samples[sample]['primerF']]
    if not single:
        primers_lgth = [samba_samples[sample]['primerF'], samba_samples[sample]['primerR']]
    cutadapt_overlap = len(min(primers_lgth, key=len)) - 1

    ## 4.3 - get the header (vars)
    if single:
        metadata_header = metadata.row_values(0, start_colx=3, end_colx=None)
    else:
        metadata_header = metadata.row_values(0, start_colx=4, end_colx=None)

    # 5 - return data
    return samba_samples, cutadapt_overlap, metadata_header

def write_qiime_input(samba_samples, metadata_header, single=False):
    manifest = open('q2_manifest.sort', 'w')
    metadata = open('q2_metadata.sort', 'w')
    # write header
    p_manifest = ['sample-id','forward-absolute-filepath','reverse-absolute-filepath']
    s_manifest = ['sample-id','absolute-filepath']
    p_metadata = ['sampleid','barcode','PrimerF','PrimerR']
    s_metadata = ['sampleid','barcode','PrimerF']
    if single:
        manifest.write('\t'.join(s_manifest) + '\n')
        metadata.write('\t'.join(s_metadata) + '\t' + '\t'.join(metadata_header) + '\n')
    else:
        manifest.write('\t'.join(p_manifest) + '\n')
        metadata.write('\t'.join(p_metadata) + '\t' + '\t'.join(metadata_header) + '\n')
    # write core file
    for sample, values in samba_samples.items():
        manifest.write(sample+'\t'+values['R1']+'\t'+values['R2']+'\n')
        metadata.write(sample+'\t'+values['barcode']+'\t'+values['primerF']+'\t'+values['primerR']+'\t'+'\t'.join(map(str,values['vars']))+'\n')
    # close
    manifest.close()
    metadata.close()

def main(args):
    # single or paired
    single=False
    if args.readtype == 'single':
        single=True

    # Manifest/metadata loading
    # If Excel input
    if args.xlsx:
        # 1 - Collect manifest and metadata infos
        print("Step 1 - load and check manifest/metadata")
        samba_samples, cutadapt_overlap, metadata_header = import_excel(args.xlsx, single=single)
    # Else, raw manifest/metadata input
    else:
        # 1.1 - Collect metadata infos and check variable names
        print("Step 1.1 - load and check q2_metadata")
        metadata, cutadapt_overlap = import_metadata(args.metadata, args.readtype)

        # 1.2 - Collect reads location from manifest
        print("Step 1.2 - load and check q2_manifest")
        manifest = import_manifest(args.manifest, args.readtype)

        # 1.3 - Check metadata vs manifest consistency
        print("Step 1.3 - check metadata vs manifest consistency")
        samba_samples = check_manifest_vs_metadata(manifest, metadata)

    # 2 - Check fastq integrity
    print("Step 2 - check fastq integrity")
    pool = mp.Pool(args.threads)
    fastq = pool.starmap(check_fastq, [(samba_samples, sample, args.readtype) for sample in samba_samples.keys()])
    pool.close()
    # Clean all this mess
    integrity_data = {}
    for result in fastq:
        for sample, vals in result.items():
            integrity_data[sample] = vals

    # 3 - Report fastq integrity before validation
    # Allow exploration by user for a better understanding
    print("Step 3 - write integrity report")
    write_report(integrity_data, args.readtype)

    # 4 - Validate manifest, metadata and fastq
    print("Step 4 - validation of data")
    integrity_validation(integrity_data, args.control, args.readtype, args.primer)

    # 5 - Write manifest and metadata or sort
    if args.xlsx:
        print("Step 5 - write manifest and metadata into Qiime2 format")
        write_qiime_input(samba_samples, metadata_header, single=single)
    else:
        print("Step 5 - Sort manifest and metadata by sample id")
        print("\tsort manifest...")
        sorted_manifest = args.manifest + ".sort"
        sort_process(args.manifest, sorted_manifest)
        print("\tsort metadata...")
        sorted_metadata = args.metadata + ".sort"
        sort_process(args.metadata, sorted_metadata)

    # 6 - Output cutadapt overlap length
    print("Step 6 - save Cutadapt overlap size")
    cutovl = open('cutovl.txt','w')
    cutovl.write(str(cutadapt_overlap))
    cutovl.close()

if __name__ == '__main__':
    args = getArgs()
    main(args)