#!/usr/bin/env python3

from loguru import logger
import rich_click as click
import json
import csv

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.version_option("1.0", prog_name="figaro_json2csv.py")
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--figaro', '-f', type=click.Path(exists=True), required=True,
              help='Output JSON from FIGARO')
def main(figaro):
    """Extract and export FIGARO results for user and nf-module"""
    # 1 - Load data
    logger.info(f'Load {figaro}')
    trim_params = load_json(figaro)
    # 2 - Export JSON to CSV
    logger.info(f'Convert FIGARO JSON to CSV')
    export_to_csv(trim_params)
    logger.success(f'Done')


def load_json(jsonfile):
    with open(jsonfile) as json_data:
        trim_params = json.load(json_data)
    return trim_params


def export_to_csv(trim_params):
    fields = ['readRetentionPercent', 'score', 'trimPosition_R1', 'trimPosition_R2', 'maxExpectedError_R1',
              'maxExpectedError_R2']

    with open('figaro.csv', 'w', newline='') as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fields)
        writer.writeheader()

        for record in trim_params:
            record['trimPosition_R1'] = record['trimPosition'][0]
            record['trimPosition_R2'] = record['trimPosition'][1]
            del record['trimPosition']
            record['maxExpectedError_R1'] = record['maxExpectedError'][0]
            record['maxExpectedError_R2'] = record['maxExpectedError'][1]
            del record['maxExpectedError']
            writer.writerow(record)


if __name__ == "__main__":
    main()
