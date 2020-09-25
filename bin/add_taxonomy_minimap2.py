#!/usr/bin/env python3

import argparse
import re
import os
import pysam

def getArgs():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-p',dest="bam",type=str,required=True,help='Path to folder with bam files from Minimap2')
    parser.add_argument('-t',dest="tax",type=str,required=True,help='Silva taxonomy database')
    parser.add_argument('-r',dest="rank",type=int,default=5,help='Minimal rank level to keep a hit as assigned [%(default)s]; 1:Kingdom, 2:Phylum, 3:Class, 4:Order, 5:Family, 6:Genus, 7:Species')
    parser.add_argument('-o',dest="out",type=str,required=True,help='Output file')

    arg = parser.parse_args()

    return arg

def store_taxonomy(taxonomy):
    taxo_r = open(taxonomy, 'r')
    taxonomy = {}
    for l in taxo_r:
        silva_seq_id, silva_tax_id = re.split(r'\t', l.rstrip('\n'))
        taxonomy[silva_seq_id] = silva_tax_id
    return taxonomy

def hits_to_taxo(taxonomy, bam_lst):
    lr_reads = {}
    for f in bam_lst:
        bam = pysam.AlignmentFile(f, "rb")
        sample = os.path.splitext(os.path.basename(f))[0]
        lr_reads[sample] = {}
        ################################################################
        # DEV: add query/subject coverage filtering (at least for 12S) #
        ################################################################
        print(sample)
        print("query_length\tasm_length\tsubject_length\tquery_coverage\tsubject_coverage")
        ################################################################
        for l in bam.fetch(until_eof=True):
            # print(l.tostring(bam))
            # Step 1 - Parse bam and collect mapped/unmapped
            ## Unmapped reads are set to no-hit / unassigned
            read = {'lr_seq_length': l.query_length }
            if l.flag == 4:
                read['lr_nm_score'] = -1
                # lr_seq_id = "no-hit"
                read['lr_tax'] = "unassigned"
            else:
                read['lr_nm_score'] = l.get_tag('NM')
                read['lr_tax'] = taxonomy[l.reference_name]
                read['lr_asm_length'] = l.query_alignment_length

                ################################################################
                # DEV: add query/subject coverage filtering (at least for 12S) #
                ################################################################
                query_length = l.infer_read_length()
                asm_length = l.query_alignment_length
                subject_length =  bam.get_reference_length(l.reference_name)
                subject_coverage = (asm_length * 100 / subject_length)
                query_coverage = (asm_length * 100 / query_length)
                printer = [str(query_length), str(asm_length), str(subject_length), str(round(query_coverage, 2)), str(round(subject_coverage, 2))]
                # print(query_length+'\t'+asm_length+'\t'+subject_length+'\t'+round(query_coverage, 2)+'\t'+round(subject_coverage, 2))
                print('\t'.join(printer))

                ################################################################

            # Step 2 - Handle new hit and multihits
            ## Add a new hit
            if not l.query_name in lr_reads[sample]:
                lr_reads[sample][l.query_name] = read
            ## Replace taxonomy if the alignment score of another hit is better
            elif lr_reads[sample][l.query_name]['lr_nm_score'] < read['lr_nm_score']:
                lr_reads[sample][l.query_name] = read
            ## In case of bame alignment score...
            elif lr_reads[sample][l.query_name]['lr_nm_score'] == read['lr_nm_score']:
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
    return lr_reads

def write_taxo(hits_to_taxo, rank, outfile):
    # sample hit count
    sample_count = ['0'] * len(hits_to_taxo)
    # Header for samples
    sample_header = []
    for s in sorted(hits_to_taxo):
        sample_header.append(s)
    # Open and start to write
    taxify = open(outfile, 'w')
    taxify.write('Read_id\t'+'\t'.join(sample_header)+'\tTaxonomy\tAssignation\tDepth\n')
    # Iterator for the position of each bample
    i = 0
    for sample in sorted(hits_to_taxo):
        # Set count to 1 for the current sample
        sample_count[i] = '1'
        for read in sorted(hits_to_taxo[sample]):
            lr_read_name = read
            lr_tax_id = hits_to_taxo[sample][lr_read_name]['lr_tax']
            lr_tax_depth = len(lr_tax_id.split(';'))
            if lr_tax_depth >= rank:
                taxify.write(lr_read_name+'\t'+'\t'.join(sample_count)+'\t'+lr_tax_id+'\tAssigned\t'+str(lr_tax_depth)+'\n')
            else:
                taxify.write(lr_read_name+'\t'+'\t'.join(sample_count)+'\t'+lr_tax_id+'\tUnassigned\t'+str(lr_tax_depth)+'\n')
        # Clean count and increment
        sample_count[i] = '0'
        i += 1

def main(args):

    # 1 - read and store taxonomy
    taxonomy = store_taxonomy(args.tax)

    # 2 - Collect all bam files
    collect_bam = [ os.path.join(args.bam, f) for f in os.listdir(args.bam) if f.endswith('.bam') ]

    # 3 - For each bam, get the taxonomy of hits
    hits_2_taxo = hits_to_taxo(taxonomy, collect_bam)

    # 4 - Filter the taxonomy by rank level and print
    # write_taxo(hits_2_taxo, args.rank, args.out)

if __name__ == '__main__':
    args = getArgs()
    main(args)
