#!/usr/bin/env python3

import xlrd
from loguru import logger
import rich_click as click


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.version_option("1.0", prog_name="excel2tsv")
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--xls', '-x', type=click.Path(exists=True), required=True,
              help='Excel file for manifest & metadata in XLS format. XLSX format is not supported.')
@click.option('--seqtype', '-s', type=click.Choice(['single', 'paired']), default='paired',
              help='Sequencing type [paired]')
def main(xls, seqtype):
    """Control and convert Excel input as metadata.tsv and manifest.tsv"""
    data_extraction = {}
    # 1 - Load data
    logger.info(f'Load {xls}')
    manifest, metadata = load_xls(xls)
    # 2 - Check file and data
    logger.info(f'Start to validate XLS')
    check_samples(manifest, metadata)
    check_seqtype(manifest, seqtype)
    check_metadata(metadata, seqtype)
    logger.success(f'Successfully validate XLS')
    # 3 - Export XLS to TSV for Qiime2
    logger.info(f'Start to export XLS to TSV')
    data_extraction = extract_manifest(manifest, seqtype, data_extraction)
    data_extraction, metadata_vars = extract_metadata(metadata, seqtype, data_extraction)
    export_to_tsv_for_qiime(data_extraction, metadata_vars, seqtype)
    logger.success(f'Done')


def load_xls(xls):
    try:
        data = xlrd.open_workbook(xls)
        manifest = data.sheet_by_name('manifest')
        metadata = data.sheet_by_name('metadata')

        logger.success(f'{xls} successfully loaded')

        return manifest, metadata

    except xlrd.XLRDError as e:
        logger.error(str(e))
        exit(1)


def check_samples(manifest, metadata):
    logger.info(f'Check for name and number of samples consistency')

    manifest_samples = set()
    metadata_samples = set()
    row_count_manifest = manifest.nrows
    row_count_metadata = metadata.nrows

    logger.info(f'Collect samples names in manifest and metadatas sheet')
    for row in range(1, row_count_manifest):
        manifest_samples.add(manifest.row_values(row, start_colx=0, end_colx=None)[0].strip(' '))
    for row in range(1, row_count_metadata):
        metadata_samples.add(metadata.row_values(row, start_colx=0, end_colx=None)[0].strip(' '))
    # Compare both to ensure samples are consistent between sheets
    if manifest_samples != metadata_samples:
        logger.error(f'Manifest and metadata did not have the same sample(s)')
        uniq_manifest = manifest_samples.difference(metadata_samples)
        uniq_metadata = metadata_samples.difference(manifest_samples)
        if len(uniq_manifest) > 0:
            logger.error(f"Sample(s) only found in manifest -> {','.join(uniq_manifest)}")
        if len(uniq_metadata) > 0:
            logger.error(f"Sample(s) only found in metadata -> {','.join(uniq_metadata)}")
        exit(1)
    else:
        logger.success(f'Manifest and metadata shared the same samples')


def check_seqtype(manifest, seqtype):
    logger.info(f'Check sequence type consistency (SE/PE)')
    col_expected = 3
    if seqtype == 'single':
        col_expected = 2
    col_count_manifest = manifest.ncols
    if col_count_manifest != col_expected:
        logger.error(
            f'Manifest did not have the required number of columns: {str(col_count_manifest)} instead of {col_expected}')
        logger.error(f'Please, configure --seqtype according to your data')
        exit(1)
    else:
        logger.success(f'Manifest is consistent with sequence type (SE/PE)')


def check_metadata(metadata, seqtype):
    logger.info(f'Check metadata formatting')
    row_count_metadata = metadata.nrows
    for row in range(1, row_count_metadata):
        tmp = metadata.row_values(row, start_colx=0, end_colx=None)
        sample = tmp[0].strip(' ')
        if seqtype == 'single':
            meta_var = [x.strip(' ') for x in map(str, tmp[3:])]
        else:
            meta_var = [x.strip(' ') for x in map(str, tmp[4:])]
        for var in meta_var:
            if var in ['NA', '', 'Nan', 'NaN', 'nan']:
                logger.error(f'{sample} has unset values (NA, NaN, empty cell...) in metadata')
                exit(1)
    logger.success(f'Metadata format is good')


def extract_manifest(manifest, seqtype, data_extraction):
    logger.info(f'Extract manifest')
    row_count_manifest = manifest.nrows
    for row in range(1, row_count_manifest):
        if seqtype == 'single':
            sample, R1 = manifest.row_values(row, start_colx=0, end_colx=None)
            data_extraction[sample] = {}
            data_extraction[sample]['R1'] = R1.strip(' ')
        else:
            sample, R1, R2 = manifest.row_values(row, start_colx=0, end_colx=None)
            data_extraction[sample] = {}
            data_extraction[sample]['R1'] = R1.strip(' ')
            data_extraction[sample]['R2'] = R2.strip(' ')

    logger.success(f'Successfully extract manifest')
    return data_extraction


def extract_metadata(metadata, seqtype, data_extraction):
    logger.info(f'Extract metadata')
    row_count_metadata = metadata.nrows
    metadata_vars = []
    for row in range(1, row_count_metadata):
        tmp = metadata.row_values(row, start_colx=0, end_colx=None)
        sample = tmp[0].strip(' ')
        data_extraction[sample]['barcode'] = tmp[1].strip(' ')
        data_extraction[sample]['primerF'] = tmp[2].strip(' ')
        if seqtype == 'single':
            metadata_vars = [x.strip(' ') for x in map(str, metadata.row_values(0, start_colx=3, end_colx=None))]
            data_extraction[sample]['vars'] = [x.strip(' ') for x in map(str, tmp[3:])]
        else:
            metadata_vars = [x.strip(' ') for x in map(str, metadata.row_values(0, start_colx=4, end_colx=None))]
            data_extraction[sample]['primerR'] = tmp[3].strip(' ')
            data_extraction[sample]['vars'] = [x.strip(' ') for x in map(str, tmp[4:])]
        for var in data_extraction[sample]['vars']:
            if var in ['NA', '', 'Nan', 'NaN', 'nan']:
                logger.error(f'{sample} has unset values (NA, NaN, empty cell...) in metadata')
                exit(1)

    logger.success(f'Successfully extract metadata')
    return data_extraction, metadata_vars


def export_to_tsv_for_qiime(data, metadata_vars, seqtype):
    logger.info(f'Start exporting XLSX to TSV for Qiime')

    manifest = open('q2_manifest.sort.tsv', 'w')
    metadata = open('q2_metadata.sort.tsv', 'w')

    p_manifest = ['sample-id', 'forward-absolute-filepath', 'reverse-absolute-filepath']
    s_manifest = ['sample-id', 'absolute-filepath']
    p_metadata = ['sampleid', 'barcode', 'PrimerF', 'PrimerR']
    s_metadata = ['sampleid', 'barcode', 'PrimerF']

    tab = '\t'
    if seqtype == 'single':
        manifest.write(f"{tab.join(s_manifest)}\n")
        metadata.write(f"{tab.join(s_metadata)}\t{tab.join(metadata_vars)}\n")
    else:
        manifest.write(f"{tab.join(p_manifest)}\n")
        metadata.write(f"{tab.join(p_metadata)}\t{tab.join(metadata_vars)}\n")
    for sample, values in data.items():
        manifest.write(f"{sample}\t{values['R1']}\t{values['R2']}\n")
        metadata.write(f"{sample}\t{values['barcode']}\t{values['primerF']}\t{values['primerR']}\t{tab.join(map(str, values['vars']))}\n")
    manifest.close()
    metadata.close()
    logger.success(f'Successfully written to TSV')


if __name__ == "__main__":
    main()
