#!/usr/bin/env python3


from __future__ import print_function
import re
import os
import pysam
from loguru import logger
import rich_click as click

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.version_option("2.0", prog_name="lr_count_table_minimap2.py")
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-b', '--bams', type=click.Path(exists=True), required=True,
              help='Path to folder with bam files from Minimap2')
@click.option('-t', '--tax', type=click.Path(exists=True), required=True, help='Silva taxonomy database')
@click.option('-r', '--rank', type=click.Choice(['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']),
              default='Family', help='Minimal rank level to keep a hit as assigned [Family]')
@click.option('-i', '--identity', type=click.FloatRange(0, 1), default=0.90, help='Minimal alignment identity [0.9]')
@click.option('-c', '--coverage', type=click.FloatRange(0, 1), default=0.90, help='Minimal query coverage [0.9]')
@click.option('-o', '--out', type=str, required=True, help='Output file')
def main(bams, tax, rank, out, identity, coverage):
    ranks = {"Kingdom": 0, "Phylum": 1, "Class": 2, "Order": 3, "Family": 4, "Genus": 5, "Species": 6}

    # 1 - read and store taxonomy
    taxonomy = store_taxonomy(tax)

    # 2 - Collect all bam files
    collect_bam = [os.path.join(bams, f) for f in os.listdir(bams) if f.endswith('.bam')]
    if len(collect_bam) == 0:
        logger.error(f'No BAM found in {bams}')
        exit(1)

    # 3 - For each bam, get the taxonomy of hits
    hits_2_taxo = hits_to_taxo(taxonomy, collect_bam, identity, coverage)

    # 4 - Filter the taxonomy by rank level and print
    write_taxo(hits_2_taxo, out, rank, ranks)

    # 5 - Debug
    # write_multihit(hits_2_taxo)


def store_taxonomy(taxonomy):
    taxo_r = open(taxonomy, 'r')
    taxonomy = {}
    for l in taxo_r:
        silva_seq_id, silva_tax_id = re.split(r'\t', l.rstrip('\n'))
        taxonomy[silva_seq_id] = silva_tax_id
    return taxonomy


def hits_to_taxo(taxonomy, bam_lst, identity, coverage):
    lr_reads = {}
    for f in bam_lst:
        bam = pysam.AlignmentFile(f, "rb")
        sample = os.path.splitext(os.path.basename(f))[0]
        lr_reads[sample] = {}
        # DEBUG
        # lr_multihit[sample] = {}
        for l in bam.fetch(until_eof=True):
            # Step 1 - Parse bam and collect mapped/unmapped
            read = {'lr_seq_length': l.query_length}
            if l.flag == 4:
                read['lr_nm_score'] = -1
                read['lr_tax'] = ''
                read['type'] = 'Unmapped'
                read['blast_like_identity'] = 'NA'
                read['query_coverage'] = 'NA'
            else:
                read['type'] = 'Assigned'
                read['lr_tax'] = taxonomy[l.reference_name]
                read['lr_nm_score'] = l.get_tag('NM')
                read['lr_asm_length'] = l.query_alignment_length
                aln_length = sum([event[1] for event in l.cigar if event[0] != 4])
                read['blast_like_identity'] = (aln_length - l.get_tag('NM')) / aln_length
                read['query_coverage'] = l.query_alignment_length / l.infer_read_length()
                read['matching_bases'] = aln_length - l.get_tag('NM')
                if read['blast_like_identity'] < identity or read['query_coverage'] < coverage:
                    read['type'] = 'Ambiguous'
            # Step 2 - Handle new hit and multihits
            # Add a new hit
            if l.query_name not in lr_reads[sample]:
                lr_reads[sample][l.query_name] = read
            # Replace hit if the alignment score of another hit is better
            elif lr_reads[sample][l.query_name]['matching_bases'] < read['matching_bases']:
                lr_reads[sample][l.query_name] = read
            # In case of same alignment score...
            elif lr_reads[sample][l.query_name]['matching_bases'] == read['matching_bases']:
                lr_reads[sample][l.query_name]['type'] = 'Multihit'
                lr_reads[sample][l.query_name]['multihit'] = []
                if len(lr_reads[sample][l.query_name]['multihit']) == 0:
                    lr_reads[sample][l.query_name]['multihit'].append(lr_reads[sample][l.query_name]['lr_tax'])
                lr_reads[sample][l.query_name]['multihit'].append(read['lr_tax'])
                # 1 - split taxo
                # c for current and n for new
                c = lr_reads[sample][l.query_name]['lr_tax'].split(';')
                n = read['lr_tax'].split(';')
                # 2 - Get tax length
                tax_depth_c = len(c)
                tax_depth_n = len(n)
                # 3 - look at the lowest/deepest rank in common
                # and store the truncated taxonomy
                for t in range(min(tax_depth_c, tax_depth_n)):
                    if c[t] != n[t]:
                        lr_reads[sample][l.query_name]['lr_tax'] = ';'.join(c[:t])
                        break
    return lr_reads


def write_taxo(hit_to_taxo, outfile, rank, ranks):
    tab = '\t'
    # sample hit count
    sample_count = ['0'] * len(hit_to_taxo)
    # Header for samples
    sample_header = []
    for s in sorted(hit_to_taxo):
        sample_header.append(s)
    # Open and start to write
    taxify = open(outfile, 'w')
    taxify.write(f'Read_id\t{tab.join(sample_header)}\tTaxonomy\tAssignation\tIdentity\tCoverage\n')
    # Iterator for the position of each sample
    i = 0
    for sample in sorted(hit_to_taxo):
        # Set count to 1 for the current sample
        sample_count[i] = '1'
        for read in sorted(hit_to_taxo[sample]):
            lr_read_name = read
            lr_tax = hit_to_taxo[sample][lr_read_name]['lr_tax']
            lr_tax_depth = len(lr_tax.split(';'))
            lr_type = hit_to_taxo[sample][lr_read_name]['type']
            lr_identity = hit_to_taxo[sample][lr_read_name]['blast_like_identity']
            lr_coverage = hit_to_taxo[sample][lr_read_name]['query_coverage']
            # Ensure the minimal rank to be 'Assigned' is fulfilled
            rank_level = ranks[rank]
            if len(lr_tax.split(';')) < rank_level and lr_type != 'Unmapped':
                lr_type = 'Ambiguous'
            # Writing
            taxify.write(f'{lr_read_name}\t{tab.join(sample_count)}\t{lr_tax}\t{lr_type}\t{lr_identity}\t{lr_coverage}\n')
        # Clean count and increment
        sample_count[i] = '0'
        i += 1
    taxify.close()


def write_multihit(hit_to_taxo):
    hits = open('multihit.txt', 'w')
    for sample in sorted(hit_to_taxo):
        for read in sorted(hit_to_taxo[sample]):
            if 'multihit' in hit_to_taxo[sample][read]:
                for hit in hit_to_taxo[sample][read]['multihit']:
                    hits.write(f'{sample}\t{read}\t{hit}\n')
    hits.close()


if __name__ == '__main__':
    main()
