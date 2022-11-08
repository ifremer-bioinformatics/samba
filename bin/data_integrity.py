#!/usr/bin/env python3

from __future__ import print_function
import re
import os
import gzip
# import subprocess
import multiprocessing as mp
from Bio import SeqIO
from loguru import logger
import rich_click as click

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.version_option("2.0", prog_name="data_integrity.py")
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--manifest', '-a', type=click.Path(exists=True), required=True, help='q2_manifest file')
@click.option('--metadata', '-e', type=click.Path(exists=True), required=True, help='q2_metadata file')
@click.option('--primer', '-p', type=click.IntRange(0, 100), default=70,
              help='Percentage of primers supposed to be found in raw reads [70]')
@click.option('--control', '-c', type=str, default='', help='Delimited list of control (comma separated)')
@click.option('--seqtype', '-s', type=click.Choice(['single', 'paired']), default='paired',
              help='Sequencing type [paired]')
@click.option('--thread', '-t', type=click.IntRange(1, max=None), default=4,
              help='Number of threads [4]')
def main(manifest, metadata, primer, control, seqtype, thread):

    extract_metadata = {}
    extract_manifest = {}

    logger.info(f"Load and check q2_metadata")
    extract_metadata = import_metadata(metadata, seqtype, extract_metadata)

    logger.info(f'Start internal checklist for q2_metadata')
    check_metadata_vars(extract_metadata)
    check_primer_sequence_uniq(extract_metadata, seqtype)

    logger.info(f"Load q2_manifest")
    extract_manifest = import_manifest(manifest, seqtype, extract_manifest)

    logger.info(f'Start internal checklist for q2_manifest')
    check_manifest_fastq_path(extract_manifest, seqtype)

    logger.info(f"Check consistency between q2_manifest and q2_metadata")
    check_manifest_vs_metadata(extract_manifest, extract_metadata)
    merged_data = merge_manifest_and_metadata(extract_manifest, extract_metadata)

    logger.info(f'Start fastq checklist from q2_manifest')
    merged_data = check_fastq_threading(merged_data, thread, seqtype)
    merged_data = convert_check_fastq_threading(merged_data)
    write_report(merged_data, seqtype)

    cutadapt_overlap = get_primer_overlap(extract_metadata, seqtype)

    logger.info('Done')


    #
    # logger.info(f"Validation of data")
    # integrity_validation(integrity_data, control, seqtype, primer)
    #
    # logger.info(f"Sort manifest and metadata by sample id")
    # logger.info(f"sort manifest")
    # sorted_manifest = manifest + ".sort"
    # sort_process(manifest, sorted_manifest)
    # logger.info(f"Sort metadata...")
    # sorted_metadata = metadata + ".sort"
    # sort_process(metadata, sorted_metadata)
    #
    # logger.info(f"Save Cutadapt overlap size")
    # cutovl = open('cutovlsize', 'w')
    # cutovl.write(str(cutadapt_overlap))
    # cutovl.close()


def import_metadata(metadata_file, seqtype, extract_metadata):
    logger.info('Import q2_metadata')
    try:
        for line in open(metadata_file, 'r'):
            if not line.startswith('sampleid'):
                sampleid = re.split(r'\t', line.rstrip('\n'))[0].strip(' ')
                extract_metadata[sampleid] = {}
                extract_metadata[sampleid]['barcode'] = re.split(r'\t', line.rstrip('\n'))[1].strip(' ')
                extract_metadata[sampleid]['primerF'] = re.split(r'\t', line.rstrip('\n'))[2].strip(' ')
                if seqtype == 'paired':
                    extract_metadata[sampleid]['primerR'] = re.split(r'\t', line.rstrip('\n'))[3].strip(' ')
                    extract_metadata[sampleid]['vars'] = [x.strip(' ') for x in map(str, re.split(r'\t', line.rstrip('\n'))[4:])]
                else:
                    extract_metadata[sampleid]['vars'] = [x.strip(' ') for x in map(str, re.split(r'\t', line.rstrip('\n'))[3:])]
        logger.success('Successfully import q2_metadata')
        return extract_metadata
    except Exception as e:
        logger.error(str(e))
        exit(1)


def check_metadata_vars(extract_metadata):
    logger.info('Check for missing or forbidden variable values in q2_metadata')
    number_var = set()
    for sample in extract_metadata:
        number_var.add(len(extract_metadata[sample]['vars']))
        for var in extract_metadata[sample]['vars']:
            if var in ['NA', '', 'Nan', 'NaN', 'nan']:
                logger.error(f'{sample} has NA/Nan/NaN/nan of empty variable value(s) in q2_metadata')
                exit(1)
    if len(number_var) != 1:
        logger.error(f'The number of variable is different between samples in q2_metadata')
        exit(1)
    logger.success('No discrepancy between variable number among samples in q2_metadata')
    logger.success('No forbidden variable values in q2_metadata')


def get_primer(extract_metadata, primer_rf):
    primer_lst = set()
    try:
        for sample in extract_metadata:
            primer_lst.add(extract_metadata[sample][primer_rf])
        return primer_lst
    except Exception as e:
        logger.error(str(e))
        exit(1)


def check_primer_sequence_uniq(extract_metadata, seqtype):
    logger.info('Check primers in q2_metadata')
    primerF = get_primer(extract_metadata, 'primerF')
    if len(primerF) != 1:
        logger.error('Multiple forward primers have been found for primerF in q2_metadata')
        exit(1)
    if seqtype == 'paired':
        primerR = get_primer(extract_metadata, 'primerR')
        if len(primerR) != 1:
            logger.error('Multiple reverse primers have been found for primerR in q2_metadata')
            exit(1)
    logger.success('No multiple primers in q2_metadata')


def get_primer_overlap(extract_metadata, seqtype):
    logger.info(f'Calculate cutadapt overlap')
    cutadapt_overlap = None
    for sample in extract_metadata:
        primers_lgth = [extract_metadata[sample]['primerF']]
        if seqtype == 'paired':
            primers_lgth = [extract_metadata[sample]['primerF'], extract_metadata[sample]['primerR']]
        cutadapt_overlap = len(min(primers_lgth, key=len)) - 1
    logger.success(f'Cutadpat overlap: {cutadapt_overlap}')
    return cutadapt_overlap


def import_manifest(manifest_file, seqtype, extract_manifest):
    logger.info('Import q2_metadata')
    try:
        for line in open(manifest_file, 'r'):
            if not line.startswith('sample-id'):
                if seqtype == 'paired':
                    try:
                        sample, R1, R2 = [x.strip(' ') for x in re.split(r'\t', line.rstrip('\n'))]
                        extract_manifest[sample] = {}
                        extract_manifest[sample]['R1'] = R1
                        extract_manifest[sample]['R2'] = R2
                    except ValueError:
                        size = len(re.split(r'\t', line.rstrip('\n')))
                        logger.error(f'q2_manifest contains {str(size)} column(s) instead of 3')
                        exit(1)
                else:
                    try:
                        sample, R1 = [x.strip(' ') for x in re.split(r'\t', line.rstrip('\n'))]
                        extract_manifest[sample] = {}
                        extract_manifest[sample]['R1'] = R1
                    except ValueError:
                        size = len(re.split(r'\t', line.rstrip('\n')))
                        logger.error(f'q2_manifest contains {str(size)} column(s) instead of 2')
                        exit(1)
        logger.success('Successfully import q2_manifest')
        return extract_manifest
    except Exception as e:
        logger.error(str(e))
        exit(1)


def check_manifest_fastq_path(extract_manifest, seqtype):
    logger.info(f'Check path to fastq')
    for sample in extract_manifest:
        check_fastq_path(extract_manifest[sample]['R1'], sample)
        if seqtype == 'paired':
            check_fastq_path(extract_manifest[sample]['R2'], sample)
    logger.success(f'Path to fastq are good')


def check_manifest_vs_metadata(manifest, metadata):
    logger.info(f'Samples names and number consistency')
    manifest_samples = set()
    metadata_samples = set()
    for k in manifest.keys():
        manifest_samples.add(k)
    for k in metadata.keys():
        metadata_samples.add(k)
    if manifest_samples != metadata_samples:
        uniq_manifest = manifest_samples.difference(metadata_samples)
        uniq_metadata = metadata_samples.difference(manifest_samples)
        logger.error(f'Manifest and metadata did not have the same sample(s)')
        if len(uniq_manifest) > 0:
            logger.error(f"Sample(s) only found in manifest -> {','.join(uniq_manifest)}")
        if len(uniq_metadata) > 0:
            logger.error(f"Sample(s) only found in metadata -> {','.join(uniq_metadata)}")
        exit(1)
    logger.success(f'Samples names and number of samples are consistent between q2_manifest and q2_metadata')


def merge_manifest_and_metadata(manifest, metadata):
    data = manifest.copy()
    for k in metadata.keys():
        for sk in metadata[k]:
            data[k][sk] = metadata[k][sk]
    return data


def check_fastq_threading(merged_data, thread, seqtype):
    logger.info("Check fastq integrity")
    pool = mp.Pool(thread)
    check_fastq_out = pool.starmap(check_fastq, [(merged_data, sample, seqtype) for sample in merged_data.keys()])
    pool.close()
    return check_fastq_out


def convert_check_fastq_threading(check_fastq_out):
    check_fastq_results = {}
    for result in check_fastq_out:
        for sample, vals in result.items():
            check_fastq_results[sample] = vals
    return check_fastq_results


def check_fastq(merged_data, sample, data_type):
    out = {}
    logger.info(f"Start analyse sample: {sample}")
    barcode = merged_data[sample]['barcode']
    reads = ['R1']
    if data_type == 'paired':
        reads.append('R2')
    for r in reads:
        R = merged_data[sample][r]
        if r == 'R1':
            primer = merged_data[sample]['primerF']
        else:
            primer = merged_data[sample]['primerR']
        instrument, index, reads_count, primer = read_fastq(R, primer)
        merged_data[sample]['reads_count_' + r] = reads_count
        merged_data[sample]['instrument_' + r] = instrument
        merged_data[sample]['nb_instrument_' + r] = len(instrument)
        merged_data[sample]['index_' + r] = index
        merged_data[sample]['primer_' + r] = primer
        merged_data[sample]['perc_primer_' + r] = round(primer * 100 / reads_count, 2)
        merged_data[sample]['nb_barcode_' + r] = len(index)
        merged_data[sample]['barcode_seq_' + r] = 'FALSE'
        if len(index) == 1:
            if index[0] == barcode:
                merged_data[sample]['barcode_seq_' + r] = 'TRUE'
    out[sample] = merged_data[sample]
    return out


def read_fastq(fastq, primer):
    primer = re.sub(r"([RYSWKMBDHVNI])", r".", primer)
    reads_count = 0
    primers_count = 0
    instrument = []
    index = []
    try:
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
    except Exception as e:
        logger.error(str(e))
        exit(1)


def check_fastq_path(path, sample):
    if not os.path.isfile(path):
        logger.error(f"Path '{path}' from q2_manifest for {sample} does not exist")
        exit(1)


def write_report(integrity_data, data_type):
    logger.info(f"Write integrity report")
    tab = '\t'
    try:
        report = open('data_integrity.txt', 'w')
        header_S = ['SampleID', 'Reads_R1', 'Barcode', 'Uniq_in_R1', 'Same_as_ref_R1', 'Uniq_sequencer_R1', 'PrimerF_in_R1',
                    'Perc_primerF_R1']
        header_P = ['SampleID', 'Reads_R1', 'Reads_R2', 'Barcode', 'Uniq_in_R1', 'Same_as_ref_R1', 'Uniq_in_R2',
                    'Same_as_ref_R2', 'Uniq_sequencer_R1', 'Uniq_sequencer_R2', 'PrimerF_in_R1', 'Perc_primerF_R1',
                    'PrimerR_in_R2', 'Perc_primerR_R2']
        if data_type == 'paired':
            report.write(f'{tab.join(header_P)}\n')
        else:
            report.write(f'{tab.join(header_S)}\n')
        for sample, val in integrity_data.items():
            if data_type == 'paired':
                report.write(
                    sample + '\t' + '{reads_count_R1}\t{reads_count_R2}\t{barcode}\t{nb_barcode_R1}\t{barcode_seq_R1}\t{nb_barcode_R2}\t{barcode_seq_R2}\t{nb_instrument_R1}\t{nb_instrument_R2}\t{primer_R1}\t{perc_primer_R1}\t{primer_R2}\t{perc_primer_R2}\n'.format(
                        **val))
            else:
                report.write(
                    sample + '\t' + '{reads_count_R1}\t{barcode}\t{nb_barcode_R1}\t{barcode_seq_R1}\t{nb_instrument_R1}\t{primer_R1}\t{perc_primer_R1}\n'.format(
                        **val))
        logger.info(f"Successfully writing report")
    except Exception as e:
        logger.error(str(e))
        exit(1)


# def integrity_validation(integrity_data, control_list, data_type, primer_threshold):
#     # collect control(s)
#     controls = [str(control) for control in control_list.split(',')]
#     # let's check the integrity
#     for sample_id, val in integrity_data.items():
#         if data_type == 'paired':
#             # 1 - check reads count R1 vs R2
#             if val['reads_count_R1'] != val['reads_count_R2']:
#                 eprint('ERROR: different number of reads between R1 and R2 for ' + sample_id)
#                 exit(1)
#             # 2 - check single sequencing instrument in reads
#             if val['nb_instrument_R1'] != 1 or val['nb_instrument_R2'] != 1:
#                 eprint('ERROR: multiple sequencing machines detected in ' + sample_id)
#                 exit(1)
#             # 3 - check primer percentage found, except for control(s)
#             if not sample_id in controls:
#                 if val['perc_primer_R1'] < primer_threshold or val['perc_primer_R2'] < primer_threshold:
#                     eprint('ERROR: ' + sample_id + " did not reach the minimum threshold for primer percentage [" + str(
#                         primer_threshold) + "%]")
#                     exit(1)
#         else:
#             # 1 - check single sequencing instrument in reads
#             if val['nb_instrument_R1'] != 1:
#                 eprint('ERROR: multiple sequencing machines detected in ' + sample_id)
#                 exit(1)
#             # 2 - check primer percentage found, except for control(s)
#             if not sample_id in controls:
#                 if val['perc_primer_R1'] < primer_threshold:
#                     eprint('ERROR: ' + sample_id + " did not reach the minimum threshold for primer percentage [" + str(
#                         primer_threshold) + "%]")
#                     exit(1)


# def sort_process(input_file, output_file):
#     cmd_sort = '(head -n 1 {in_file} && tail -n +2 {in_file} | sort) > {out_file}'.format(in_file=input_file,
#                                                                                           out_file=output_file)
#     p = subprocess.Popen(cmd_sort, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
#     stdout, stderr = p.communicate()
#     if p.returncode != 0:
#         eprint('ERROR: sort metadata/manifest failed')
#         raise Exception(stderr)
#         exit(1)


# def write_qiime_input(samba_samples, metadata_header, single=False):
#     manifest = open('q2_manifest.sort', 'w')
#     metadata = open('q2_metadata.sort', 'w')
#     # write header
#     p_manifest = ['sample-id', 'forward-absolute-filepath', 'reverse-absolute-filepath']
#     s_manifest = ['sample-id', 'absolute-filepath']
#     p_metadata = ['sampleid', 'barcode', 'PrimerF', 'PrimerR']
#     s_metadata = ['sampleid', 'barcode', 'PrimerF']
#     if single:
#         manifest.write('\t'.join(s_manifest) + '\n')
#         metadata.write('\t'.join(s_metadata) + '\t' + '\t'.join(metadata_header) + '\n')
#     else:
#         manifest.write('\t'.join(p_manifest) + '\n')
#         metadata.write('\t'.join(p_metadata) + '\t' + '\t'.join(metadata_header) + '\n')
#     # write core file
#     for sample, values in samba_samples.items():
#         manifest.write(sample + '\t' + values['R1'] + '\t' + values['R2'] + '\n')
#         metadata.write(
#             sample + '\t' + values['barcode'] + '\t' + values['primerF'] + '\t' + values['primerR'] + '\t' + '\t'.join(
#                 map(str, values['vars'])) + '\n')
#     # close
#     manifest.close()
#     metadata.close()


if __name__ == '__main__':
    main()
