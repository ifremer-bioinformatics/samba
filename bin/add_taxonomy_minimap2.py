#!/usr/bin/env python3

import argparse
import re
import os

def getArgs():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-p',dest="sam",type=str,required=True,help='Path to folder with SAM files from Minimap2')
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

def hits_to_taxo(taxonomy, sam_lst):
    # TODO: gestion des headers et du parsing en utilisant le q2_manifest?
    reads = {}
    for f in sam_lst:
        sam = open(f, 'r')
        sample = os.path.splitext(os.path.basename(f))[0]
        reads[sample] = {}
        for l in sam:
            # Step 1 - Parse sam and collect mapped/unmapped
            longreads_seq_id = l.split()[0]
            sam_flag = l.split()[1]
            silva_seq_id = l.split()[2]

            if not sam_flag == "4":
                align_score = int(l.split()[13].split(':')[2])
                silva_tax_id = taxonomy[silva_seq_id]
            else:
                align_score = -1
                silva_seq_id = "no-hit"
                silva_tax_id = "unassigned"

            # Step 2 - Handle new hit and multihits
            ## Dictionnary structure
            """
            reads: {
                sample: {
                    nanopore_read_id: [silva_tax_id, align_score]
                    ...
                }
            }
            """
            ## Add a new hit
            if not longreads_seq_id in reads[sample]:
                reads[sample][longreads_seq_id] = [silva_tax_id, align_score]

            ## Replace taxonomy if the alignment score of another hit is better
            elif reads[sample][longreads_seq_id][1] < align_score:
                reads[sample][longreads_seq_id] = [silva_tax_id, align_score]

            ## In case of same alignment score...
            elif reads[sample][longreads_seq_id][1] == align_score:
                ### 1 - split taxo
                # c for current and n for new
                c =  reads[sample][longreads_seq_id][0].split(';')
                n =  silva_tax_id.split(';')
                ### 2 - Get tax length
                tax_depth_c = len(c)
                tax_depth_n = len(n)
                ### 3 - look at the highest rank in common
                ### and store the truncated taxonomy
                for l in range(min(tax_depth_c, tax_depth_n)):
                    if c[l] != n[l]:
                        reads[sample][longreads_seq_id] = [';'.join(c[:l]), align_score]
                        break
    return reads

def write_taxo(hits_to_taxo, rank, outfile):
    # Sample hit count
    sample_count = ['0'] * len(hits_to_taxo)
    # Header for samples
    sample_header = []
    for s in sorted(hits_to_taxo):
        sample_header.append(s)
    # Open and start to write
    taxify = open(outfile, 'w')
    taxify.write('Read_id\t'+'\t'.join(sample_header)+'\tTaxonomy\tAssignation\tDepth\n')
    # Iterator for the position of each sample
    i = 0
    for sample in sorted(hits_to_taxo):
        # Set count to 1 for the current sample
        sample_count[i] = '1'
        for read in sorted(hits_to_taxo[sample]):
            longreads_seq_id = read
            silva_tax_id = hits_to_taxo[sample][longreads_seq_id][0]
            tax_depth = len(silva_tax_id.split(';'))
            if tax_depth >= rank:
                taxify.write(longreads_seq_id+'\t'+'\t'.join(sample_count)+'\t'+silva_tax_id+'\tAssigned\t'+str(tax_depth)+'\n')
            else:
                taxify.write(longreads_seq_id+'\t'+'\t'.join(sample_count)+'\t'+silva_tax_id+'\tUnassigned\t'+str(tax_depth)+'\n')
        # Clean count and increment
        sample_count[i] = '0'
        i += 1

def main(args):

    # 1 - read and store taxonomy
    taxonomy = store_taxonomy(args.tax)

    # 2 - Collect all sam files
    collect_sam = [ os.path.join(args.sam, f) for f in os.listdir(args.sam) if f.endswith('.sam') ]

    # 3 - For each sam, get the taxonomy of hits
    hits_2_taxo = hits_to_taxo(taxonomy, collect_sam)

    # 4 - Filter the taxonomy by rank level and print
    write_taxo(hits_2_taxo, args.rank, args.out)

if __name__ == '__main__':
    args = getArgs()
    main(args)
