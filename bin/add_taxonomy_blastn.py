#!/usr/bin/env python3

from __future__ import print_function
import os
import sys
import gzip
import taxopy
import argparse
from Bio import SeqIO

# For errors / warnings
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def getArgs():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-b',dest="blast",type=str,required=True,help='Blast result')
    parser.add_argument('-f',dest="fastq",type=str,required=True,help='Fastq file')
    parser.add_argument('-t',dest="taxdb",type=str,required=True,help='Path to nodes.dmp and names.dmp')
    parser.add_argument('-o',dest="output",type=str,required=True,help='Output file name')

    arg = parser.parse_args()

    return arg

def load_dmp(tax_dmp_path):
    nodes = os.path.join(tax_dmp_path, "nodes.dmp")
    names = os.path.join(tax_dmp_path, "names.dmp")
    dmp = taxopy.TaxDb(nodes_dmp=nodes, names_dmp=names, keep_files=True)
    return dmp

def parse_blast(blast_best_hit, dmp):
    '''
    outfmt="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms"
    sort -k1,1 -k12,12nr -k11,11n blast_output | sort -u -k1,1 --merge > blast_output_besthit
    '''
    taxify = {}
    for line in open(blast_best_hit, 'r'):
        query = line.split()[0]
        tax_id = line.split()[12].split(';')[0]
        try:
            lineage = taxopy.Taxon(tax_id, dmp)
            lineage_lst = lineage.name_lineage
            lineage_lst.reverse()
            taxify[query] = lineage_lst
        except taxopy.exceptions.TaxidError:
            eprint('WARNING: read '+ query + ' with taxid ' + tax_id + ' which is not present in database. Please, index your BlastDB with the same version of NCBI Taxonomy.')
    return taxify

def write_taxo(fastq, blast_taxify, output):
    blast_tax = open(output, 'w')
    with gzip.open(fastq, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            seq_id = record.description.split()[0]
            if seq_id in blast_taxify:
                blast_tax.write(seq_id + '\t' + ';'.join(blast_taxify[seq_id]) + '\n')
            else:
                blast_tax.write(seq_id + '\tNoHit\n')

def main(args):

    # 1 - load nodes and names dmp
    dmp = load_dmp(args.taxdb)

    # 2 - parse blast result
    blast_taxify = parse_blast(args.blast, dmp)

    # 3 - write the results
    write_taxo(args.fastq, blast_taxify, args.output)

if __name__ == '__main__':
    args = getArgs()
    main(args)