#!/usr/bin/env python3

import argparse
import json
import re
import os
from jinja2 import Template
from collections import defaultdict

def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-t',dest="template",type=str,required=True,help='HTML template for Jinja2')
    parser.add_argument('-c', dest="config", type=str, required=True, help='config file - JSON')
    parser.add_argument('-p', dest="path", type=str, required=True, help='PATH/TO/00_report/')

    args = parser.parse_args()

    return args

def load_json(jsonfile):
    nxt_params = {}
    with open(jsonfile) as json_data:
        nxt_params = json.load(json_data)
    return nxt_params

def count_seq_fasta(fasta):
    r_fasta = open(fasta, 'r')
    count = 0
    for l in r_fasta:
        if l.startswith('>'):
            count += 1
    return count

def count_lines_file(file):
    with open(file) as f:
        count = sum(1 for line in f)
    return count

def collect_from_folder(path):
    files = os.listdir(path)
    return files

def walk(directory):
    for root, dirs, files in os.walk(directory):
        yield root, dirs, files

def write_html(template, params, results):
    template = open(template, 'r').read()
    report = open('SAMBA_report.html', 'w')
    # Template rendering
    tm = Template(template)
    msg = tm.render(d=params, r=results)
    report.write(msg)
    # Close file
    report.close()

def get_true_from_csv(csv):
    csv_r = open(csv, 'r')
    valid_asv = []
    for l in csv_r:
        try:
            asv, value, valid = l.split()
            if valid == 'True':
                valid_asv.append(asv)
        # Header in R format is painful...
        except:
            pass
    return valid_asv

def main(args):
    # Path and files to render
    # TODO: put that into a JSON file
    structure = {
        'dada2': {
            'fasta': 'dada2_output/sequences.fasta'
        },
        'dbotu3': {
            'fasta': 'dbotu3_output/sequences.fasta'
        },
        'tax_filtering': {
            'fasta': 'tax_filtering/tax_filtered_output/sequences.fasta'
        },
        'microdecon': {
            'fasta': 'microDecon/decontaminated_ASV.fasta',
            'removed': 'microDecon/ASV_removed.txt'
        },
        'picrust': {
            'folder': 'picrust2_output/'
        },
        'ancom': {
            'folder': 'ancom_output/'
        },
        'alpha_div': {
            'folder': 'R/FIGURES/alpha_diversity/'
        },
        'beta_div': {
            'folder': 'R/FIGURES/'
        },
        'descriptive_comparison': {
            'folder': 'R/FIGURES/descriptive_comparison/'
        }

    }

    # -------------------------------------------------------------
    # -- Main -----------------------------------------------------
    # -------------------------------------------------------------
    # Store results into dict
    results = defaultdict(dict)

    # _____________________________________________________________
    # Step 1 - load the json file with all the Nextflow parameters
    nxt_params = load_json(args.config)

    # Only keep thename instead of full path
    nxt_params['metadata'] = os.path.basename(nxt_params['metadata'])
    nxt_params['manifest'] = os.path.basename(nxt_params['manifest'])

    # Store methods, indices and stats
    nxt_params['methods'] = ['rarefied','DESeq2','CSS']
    nxt_params['indices'] = ['NMDS', 'PCoA', 'hclustering']
    nxt_params['stats'] = ['bray', 'wunifrac', 'jaccard', 'unifrac']

    # _____________________________________________________________
    # Step 2 - Collect info and files for each process
    # _____________________________________________________________
    # Step 2.1 - Short reads
    # _____________________________________________________________
    if nxt_params['steps']['longreads'] == 'false':

        # Only keep the DB name instead of full path
        nxt_params['taxonomy']['database'] = os.path.basename(nxt_params['taxonomy']['database'])

        ## DADA 2
        dada2_fasta = os.path.join(args.path, structure['dada2']['fasta'])
        results['dada2']['asv_count'] = count_seq_fasta(dada2_fasta)

        # _____________________________________________________________
        ## DBOTU 3
        if nxt_params['steps']['dbotu3_enable'] == 'true':
            dbotu3_fasta = os.path.join(args.path,structure['dbotu3']['fasta'])
            results['dbotu3']['asv_count'] = count_seq_fasta(dbotu3_fasta)
            results['dbotu3']['clustering'] = results['dada2']['asv_count'] - results['dbotu3']['asv_count']
            results['dbotu3']['clustering_perc'] = round(100 - (results['dbotu3']['asv_count'] * 100.00 / results['dada2']['asv_count']), 2)

        # _____________________________________________________________
        ## Tax filtering
        if nxt_params['steps']['filtering_tax_enable'] == 'true':
            tax_filtering_fasta = os.path.join(args.path,structure['tax_filtering']['fasta'])
            results['tax_filtering']['asv_count'] = count_seq_fasta(tax_filtering_fasta)
            if nxt_params['steps']['dbotu3_enable'] == 'true':
                results['tax_filtering']['asv_rm'] = results['dbotu3']['asv_count'] - results['tax_filtering']['asv_count']
            else:
                results['tax_filtering']['asv_rm'] = results['dada2']['asv_count'] - results['tax_filtering']['asv_count']

        # _____________________________________________________________
        ## MicroDecon
        if nxt_params['steps']['microDecon_enable'] == 'true':
            microdecon_fasta = os.path.join(args.path,structure['microdecon']['fasta'])
            microdecon_removed = os.path.join(args.path,structure['microdecon']['removed'])

            results['microdecon']['asv_count'] = count_seq_fasta(microdecon_fasta)
            if nxt_params['steps']['filtering_tax_enable'] == 'true':
                results['microdecon']['asv_rm'] = results['tax_filtering']['asv_count'] - results['microdecon']['asv_count']
            elif nxt_params['steps']['dbotu3_enable'] == 'true':
                results['microdecon']['asv_rm'] = results['dbotu3']['asv_count'] - results['microdecon']['asv_count']
            else:
                results['microdecon']['asv_rm'] = results['dada2']['asv_count'] - results['microdecon']['asv_count']
            results['microdecon']['asv_rm_ctrl'] = count_lines_file(microdecon_removed) - 1

        # _____________________________________________________________
        ## piCrust
        if nxt_params['steps']['picrust2_enable'] == 'true':
            picrust_out_dir = os.path.join(args.path,structure['picrust']['folder'])
            picrust_listdir = collect_from_folder(picrust_out_dir)

            ### Collect png for each category
            results['picrust'] = defaultdict(dict)
            results['picrust']['var_tested'] = []
            for file in picrust_listdir:
                if file.startswith('EC_') and file.endswith('.png'):
                    var_name = os.path.splitext(file)[0].replace('EC_functional_predictions_NMDS_','')
                    results['picrust']['var_tested'].append(var_name)

        # _____________________________________________________________
        ## ancom
        if nxt_params['steps']['ancom_enable'] == 'true':
            ancom_out_dir = os.path.join(args.path,structure['ancom']['folder'])
            ancom_listdir = collect_from_folder(ancom_out_dir)

        ### Collect variable names
        if nxt_params['steps']['ancom_enable'] == 'true':
            results['ancom']['var'] = []
            for f in ancom_listdir:
                if not f.endswith('_family') and not f.endswith('_genus'):
                    var = re.sub(r'export_ancom_', '', f)
                    results['ancom']['var'].append(var)

        ### Collect ASV folder and DE
        if nxt_params['steps']['ancom_enable'] == 'true':
            for var in results['ancom']['var']:
                # Folder path
                var_path_base = os.path.join(ancom_out_dir, ('export_ancom_'+ var))
                var_path_family = os.path.join(ancom_out_dir, ('export_ancom_' + var + '_family'))
                var_path_genus = os.path.join(ancom_out_dir, ('export_ancom_' + var + '_genus'))
                # ASV DE
                results['ancom'][var] = {}
                results['ancom'][var]['base'] = len(get_true_from_csv(os.path.join(var_path_base, 'ancom.tsv')))
                results['ancom'][var]['family'] = len(get_true_from_csv(os.path.join(var_path_family, 'ancom.tsv')))
                results['ancom'][var]['genus'] = len(get_true_from_csv(os.path.join(var_path_genus, 'ancom.tsv')))

    # _____________________________________________________________
    # Step 2.2 - Long reads
    # _____________________________________________________________
    else:
        # _____________________________________________________________
        ## long reads databases and stats
        nxt_params["lr"]["lr_tax_fna"] = os.path.basename(nxt_params["lr"]["lr_tax_fna"])
        nxt_params["lr"]["lr_tax_flat"] = os.path.basename(nxt_params["lr"]["lr_tax_flat"])
        nxt_params['stats'] = ['bray', 'jaccard']
        if nxt_params["lr"]["lr_type"] == 'map-ont':
            nxt_params["lr"]["lr_type"] = 'Nanopore'
        else:
            nxt_params["lr"]["lr_type"] = 'PacBio'
        ranks = {'1': 'Kingdom', '2': 'Phylum', '3': 'Class', '4': 'Order', '5': 'Family', '6': 'Genus', '7': 'Species'}
        nxt_params["lr"]["lr_rank"] = ranks[nxt_params["lr"]["lr_rank"]]

    # _____________________________________________________________
    # Step 2.3 - Short and long reads
    # _____________________________________________________________
    ## Alpha diversity
    if nxt_params['steps']['stats_alpha_enable'] == 'true':
        alpha_div_out_dir = os.path.join(args.path,structure['alpha_div']['folder'])

        ### Collect all variables tested, using barplot folder
        alpha_div_barplot_out_dir = os.path.join(alpha_div_out_dir, 'diversity_barplots')
        results['alpha_div']['var_tested'] = collect_from_folder(alpha_div_barplot_out_dir)

    # _____________________________________________________________
    ## Beta diversity
    if nxt_params['steps']['stats_beta_enable'] == 'true':
        beta_div_out_dir = os.path.join(args.path,structure['beta_div']['folder'])

        # Collect all variables tested, using non_normalized folder
        non_normalized_out_dir_tested_var = os.path.join(beta_div_out_dir, 'beta_diversity_non_normalized', 'NMDS')
        results['beta_div']['var_tested'] = []
        for f in collect_from_folder(non_normalized_out_dir_tested_var):
            if f.endswith('_bray.png'):
                f = f.replace('NMDS_','').replace('_bray.png','')
                results['beta_div']['var_tested'].append(f)

    # _____________________________________________________________
    ## Descriptive comparison
    if nxt_params['steps']['stats_desc_comp_enable'] == 'true':
        descriptive_comparison_out_dir = os.path.join(args.path,structure['descriptive_comparison']['folder'])
        descriptive_comparison_listdir = collect_from_folder(descriptive_comparison_out_dir)

        # Collect files
        results['descriptive_comparison']['var'] = {}
        results['descriptive_comparison']['var_tested'] = []
        for f in descriptive_comparison_listdir:
            if f.endswith('.png'):
                var_name = f.replace('upset_plot_', '').replace('.png', '')
                results['descriptive_comparison']['var_tested'].append(var_name)

    # -------------------------------------------------------------
    # Step 3 - write the final html report
    # _____________________________________________________________
    write_html(args.template, nxt_params, results)

if __name__ == '__main__':
    args = get_args()
    main(args)
