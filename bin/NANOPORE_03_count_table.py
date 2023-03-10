#!/usr/bin/env python3


from __future__ import print_function
import re
import os
import sys
import pysam
import argparse

def getArgs():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-p',dest="bam",type=str,required=True,help='Path to folder with bam files from Minimap2')
    parser.add_argument('-t',dest="tax",type=str,required=True,help='Silva taxonomy database')
    parser.add_argument('-r',dest="rank",type=int,default=5,help='Minimal rank level to keep a hit as assigned [%(default)s]; 1:Kingdom, 2:Phylum, 3:Class, 4:Order, 5:Family, 6:Genus, 7:Species')
    parser.add_argument('-o',dest="out",type=str,required=True,help='Output file')

    arg = parser.parse_args()

    return arg

# For errors / warnings
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def store_taxonomy(taxonomy):
    taxo_r = open(taxonomy, 'r')
    taxonomy = {}
    for l in taxo_r:
        silva_seq_id, silva_tax_id = re.split(r'\t', l.rstrip('\n'))
        taxonomy[silva_seq_id] = silva_tax_id
    return taxonomy

def hits_to_taxo(taxonomy, bam_lst):
    lr_reads = {}
    lr_ambiguous_reads = {}
    for f in bam_lst:
        bam = pysam.AlignmentFile(f, "rb")
        sample = os.path.splitext(os.path.basename(f))[0]
        lr_reads[sample] = {}
        ################################################################
        # DEV: add query/subject coverage filtering (at least for 12S) #
        ################################################################
        # print(sample)
        # print("query_length\tasm_length\tsubject_length\tquery_coverage\tsubject_coverage")
        ################################################################
        for l in bam.fetch(until_eof=True):
            # print(l.tostring(bam))
            # Step 1 - Parse bam and collect mapped/unmapped
            ## Unmapped reads are set to no-hit / unassigned
            read = {'lr_seq_length': l.query_length }
            if l.flag == 4:
                read['lr_nm_score'] = -1
                # lr_seq_id = "no-hit"
                read['lr_tax'] = "Unmapped"
            else:
                read['lr_nm_score'] = l.get_tag('NM')
                read['lr_tax'] = taxonomy[l.reference_name]
                read['lr_asm_length'] = l.query_alignment_length

                # compute aln length from cigar line
                aln_length = sum([event[1] for event in l.cigar if event[0] != 4])

                # Possible metric to compute from the alignement length:

                # read['blast_like_identity'] = (aln_length - l.get_tag('NM') )/aln_length
                # read['gap_compressed_identity'] = 1 - l.get_tag('de')
                # read['query_coverage'] = l.query_alignment_length / l.infer_read_length()

                # select matching_bases because it give an idea of length and identity of the alignment.
                read['matching_bases'] = aln_length - l.get_tag('NM')

                ################################################################
                # DEV: add query/subject coverage filtering (at least for 12S) #
                ################################################################
                # query_length = l.infer_read_length()
                # asm_length = l.query_alignment_length
                # subject_length =  bam.get_reference_length(l.reference_name)
                # subject_coverage = (asm_length * 100 / subject_length)
                # query_coverage = (asm_length * 100 / query_length)
                # printer = [str(query_length), str(asm_length), str(subject_length), str(round(query_coverage, 2)), str(round(subject_coverage, 2))]
                # print('\t'.join(printer))
                ################################################################

            # Step 2 - Handle new hit and multihits
            ## Add a new hit
            if not l.query_name in lr_reads[sample]:
                lr_reads[sample][l.query_name] = read
            ## Replace taxonomy if the alignment score of another hit is better
            elif lr_reads[sample][l.query_name]['matching_bases'] < read['matching_bases']:
                lr_reads[sample][l.query_name] = read
            ## In case of same alignment score...
            elif lr_reads[sample][l.query_name]['matching_bases'] == read['matching_bases']:
                ### 1 - split taxo
                # c for current and n for new
                c =  lr_reads[sample][l.query_name]['lr_tax'].split(';')
                n =  read['lr_tax'].split(';')
                ### 2 - Get tax length
                tax_depth_c = len(c)
                tax_depth_n = len(n)
                ### 3 - look at the lowest/deepest rank in common
                ### and store the truncated taxonomy
                for t in range(min(tax_depth_c, tax_depth_n)):
                    if c[t] != n[t]:
                        # lr_reads[sample][l.query_name] = read
                        lr_reads[sample][l.query_name]['lr_tax'] = ';'.join(c[:t])
                        break
                ### 3 - Store the taxonomy options
                # if lr_ambiguous_reads:
                # else:

    return lr_reads

def write_taxo(hits_to_taxo, rank, outfile, ranks):
    # sample hit count
    sample_count = ['0'] * len(hits_to_taxo)
    # Header for samples
    sample_header = []
    for s in sorted(hits_to_taxo):
        sample_header.append(s)
    # Open and start to write
    taxify = open(outfile, 'w')
    taxify.write('Read_id\t'+'\t'.join(sample_header)+'\tTaxonomy\tAssignation\n')
    # Iterator for the position of each bample
    i = 0
    for sample in sorted(hits_to_taxo):
        # Set count to 1 for the current sample
        sample_count[i] = '1'
        for read in sorted(hits_to_taxo[sample]):
            lr_read_name = read
            lr_tax = hits_to_taxo[sample][lr_read_name]['lr_tax']
            lr_tax_depth = len(lr_tax.split(';'))
            # Add ambiguous field(s)
            if lr_tax_depth != 7:
                lr_tax = ambiguous_taxa(lr_tax, lr_tax_depth, ranks)
            # Writing
            if lr_tax_depth >= rank:
                taxify.write(lr_read_name+'\t'+'\t'.join(sample_count)+'\t'+lr_tax+'\tAssigned\n')
            elif lr_tax == "Unmapped":
                taxify.write(lr_read_name+'\t'+'\t'.join(sample_count)+'\t'+lr_tax+'\tUnmapped\n')
            else:
                taxify.write(lr_read_name+'\t'+'\t'.join(sample_count)+'\t'+lr_tax+'\tAmbiguous\n')

        # Clean count and increment
        sample_count[i] = '0'
        i += 1

def ambiguous_taxa(lr_tax, lr_tax_depth, ranks):
    ranks_depth = len(ranks)
    for i in range(lr_tax_depth, ranks_depth):
        lr_tax = lr_tax + ";Ambiguous_" + ranks[i]
    return lr_tax

# def write_ambiguous_reads():


def main(args):

    ranks = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

    # 1 - read and store taxonomy
    taxonomy = store_taxonomy(args.tax)

    # 2 - Collect all bam files
    collect_bam = [ os.path.join(args.bam, f) for f in os.listdir(args.bam) if f.endswith('.bam') ]
    if len(collect_bam) == 0:
        eprint('ERROR: no BAM found')
        exit(1)

    # 3 - For each bam, get the taxonomy of hits
    hits_2_taxo = hits_to_taxo(taxonomy, collect_bam)

    # 4 - Filter the taxonomy by rank level and print
    write_taxo(hits_2_taxo, args.rank, args.out, ranks)

if __name__ == '__main__':
    args = getArgs()
    main(args)
