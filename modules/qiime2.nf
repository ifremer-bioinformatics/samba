process q2_import_data {

    label 'qiime2_env'

    publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: '03_import_output'
    publishDir "${params.outdir}/${params.import_step}", mode: 'copy', pattern: 'data.qz*'
    publishDir "${params.outdir}/${params.report_dirname}/98_version", mode: 'copy', pattern: 'v_*.txt'
    publishDir "${params.outdir}/${params.report_dirname}/99_completecmd", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_import -> "03_${task.process}_complete.sh" }

    input:
        path(final_manifest)

    output:
        path('data.qza'), emit: imported_data
        path('data.qzv')
        path('03_import_output')
        path('completecmd')
        path('v_qiime2.txt')

    script:
    """
    03_q2_import.sh ${params.singleEnd} ${final_manifest} data.qza data.qzv 03_import_output completecmd &> q2_import.log 2>&1
    qiime --version|grep 'q2cli'|cut -d' ' -f3 > v_qiime2.txt
    """

}

process q2_cutadapt {

    label 'qiime2_env'
    label 'multithreads'

    publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: '04_cutadapt_output'
    publishDir "${params.outdir}/${params.cutadapt_step}", mode: 'copy', pattern: 'trimmed_data.qz*'
    publishDir "${params.outdir}/${params.report_dirname}/98_version", mode: 'copy', pattern : 'v_cutadapt.txt'
    publishDir "${params.outdir}/${params.report_dirname}/99_completecmd", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_cutadapt -> "cmd/04_${task.process}_complete.sh" }

    input:
        path(imported_data)

    output:
        path('trimmed_data.qza'), emit: trimmed_data
        path('trimmed_data.qzv')
        path('04_cutadapt_output')
        path('v_cutadapt.txt')
        path('completecmd')

    script:
    """
    04_q2_cutadapt.sh ${params.singleEnd} ${task.cpus} ${imported_data} ${params.primerF} ${params.primerR} ${params.errorRate} trimmed_data.qza trimmed_data.qzv 04_cutadapt_output completecmd &> q2_cutadapt.log 2>&1
    cutadapt --version > v_cutadapt.txt
    """

}

process q2_dada2 {

    label 'qiime2_env'

    publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: '06_DADA2_output'
    publishDir "${params.outdir}/${params.dada2_step}", mode: 'copy', pattern: 'DADA2_*'
    publishDir "${params.outdir}/${params.report_dirname}/98_version", mode: 'copy', pattern : 'v_dada2.txt'
    publishDir "${params.outdir}/${params.report_dirname}/99_completecmd", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_dada2 -> "cmd/06_${task.process}_complete.sh" }

    input:
        path(dada2_input_data)
        path(metadata)
        val(dada2_params)

    output:
        path('DADA2_rep_seqs.qza'), emit: dada2_rep_seqs
        path('DADA2_rep_seqs.qzv')
        path('DADA2_table.qza'), emit: dada2_table
        path('DADA2_table.qzv')
        path('DADA2_process_stats.qz*')
        path('06_DADA2_output'), emit: dada2_outdir
        path('06_DADA2_output/dada2_asv_table.tsv'), emit: dada2_table_tsv
        path('06_DADA2_output/sequences.fasta'), emit: dada2_asv_seqs_fasta
        path('sequences_with_abundance.fasta'), emit: dada2_asv_seqs_fasta_abundance
        path('v_dada2.txt')
        path('completecmd')

    script:
    """
    06_q2_dada2.sh ${params.singleEnd} ${dada2_input_data} ${params.figaro_enable} ${params.raw_read_length} ${params.primerF} ${params.primerR} ${params.FtrimLeft} ${params.RtrimLeft} ${dada2_params[0]} ${dada2_params[1]} ${params.truncQ} ${dada2_params[2]} ${dada2_params[3]} ${params.n_read_learn} ${params.pooling_method} ${params.chimeras_method} ${task.cpus} DADA2_rep_seqs.qza DADA2_table.qza DADA2_process_stats.qza DADA2_process_stats.qzv DADA2_table.qzv ${metadata} DADA2_rep_seqs.qzv 06_DADA2_output completecmd &> q2_dada2.log 2>&1
    echo '1.26.0' > v_dada2.txt
    """

}

process q2_dbOTU3 {

    label 'qiime2_env'

    publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: '07_dbOTU3_output'
    publishDir "${params.outdir}/${params.dbotu3_step}", mode: 'copy', pattern: 'dbOTU3*'
    publishDir "${params.outdir}/${params.report_dirname}/98_version", mode: 'copy', pattern: 'v_dbotu3.txt'
    publishDir "${params.outdir}/${params.report_dirname}/99_completecmd", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_dbotu3 -> "cmd/07_${task.process}_complete.sh" }

    input:
        path(asv_table)
        path(asv_seqs)
        path(metadata)

    output:
        path('dbOTU3_details.txt')
        path('dbOTU3_seqs.qza'), emit: dbotu3_seqs
        path('dbOTU3_seqs.qzv')
        path('dbOTU3_table.qza'), emit: dbotu3_table
        path('dbOTU3_table.qzv')
        path('07_dbOTU3_output'), emit: dbotu3_outdir
        path('07_dbOTU3_output/sequences.fasta'), emit: dbotu3_asv_seqs_fasta
        path('v_dbotu3.txt')
        path('completecmd')

    script:
    """
    07a_q2_dbotu3.sh ${asv_table} ${asv_seqs} ${params.gen_crit} ${params.abund_crit} ${params.pval_crit} dbOTU3_seqs.qza dbOTU3_table.qza dbOTU3_details.txt dbOTU3_table.qzv ${metadata} dbOTU3_seqs.qzv 07_dbOTU3_output completecmd &> q2_dbotu3.log 2>&1
    echo "1.5.3" > v_dbotu3.txt
    """

}

process q2_assign_taxo {

    label 'qiime2_env'
    label 'highRAM'

    publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: '08_taxonomic_assignation_output'
    publishDir "${params.outdir}/${params.assign_taxo_step}", mode: 'copy', pattern: '*.qz*'
    publishDir "${params.outdir}/${params.assign_taxo_step}", mode: 'copy', pattern: 'ASV_tax_table.biom'
    publishDir "${params.outdir}/${params.report_dirname}/99_completecmd", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_assign_taxo -> "cmd/08_${task.process}_complete.sh" }

    input:
        path(asv_seqs)
        path(asv_outdir)

    output:
        path('taxonomy.qza'), emit: taxonomy_assigned
        path('taxonomy.qzv')
        path('taxonomy_assigned.tsv'), emit: taxonomy_tsv
        path('08_taxonomic_assignation_output')
        path('ASV_tax_table.biom'), emit: asv_tax_table_biom
        path('ASV_tax_table.tsv'), emit: asv_tax_table_tsv
        path('completecmd')

    script:
    """
    08_q2_assign_taxo.sh qiime2_tmpdir ${task.cpus} ${params.confidence} ${params.database} ${asv_seqs} taxonomy.qza taxonomy.qzv 08_taxonomic_assignation_output taxonomy_assigned.tsv ${asv_outdir} ASV_tax_table.biom ASV_tax_table.tsv completecmd &> q2_assign_taxo.log 2>&1
    """

}

process q2_filter_table_by_tax {

    label 'qiime2_env'

    publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: '09_filter_table_by_tax_output'
    publishDir "${params.outdir}/${params.filter_table_by_tax_step}", mode: 'copy', pattern: '*.qz*'
    publishDir "${params.outdir}/${params.filter_table_by_tax_step}", mode: 'copy', pattern: 'asv_table_tax_filtered.biom'
    publishDir "${params.outdir}/${params.report_dirname}/99_completecmd", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_filtering_tax -> "cmd/09_${task.process}_complete.sh" }

    input:
        path(asv_table)
        path(asv_seqs)
        path(tax_qza)
        path(tax_tsv)
        path(metadata)

    output:
        path('asv_table_tax_filtered.qza'), emit: asv_table_tax_filtered_qza
        path('asv_table_tax_filtered.qzv')
        path('asv_seqs_tax_filtered.qza'), emit: asv_seqs_tax_filtered_qza
        path('asv_seqs_tax_filtered.qzv')
        path('09_filter_table_by_tax_output'), emit: filter_table_by_tax_output
        path('09_filter_table_by_tax_output/sequences.fasta'), emit: filter_table_by_tax_seqs_fasta
        path('asv_table_tax_filtered.biom'), emit: asv_table_tax_filtered_biom
        path('asv_table_tax_filtered.tsv'), emit: asv_table_tax_filtered_tsv
        path('completecmd')

    script:
    """
    09_q2_filter_table_by_tax.sh ${params.filtering_type} ${asv_table} ${tax_qza} ${params.tax_to_filter} asv_table_tax_filtered.qza ${asv_seqs} asv_seqs_tax_filtered.qza asv_table_tax_filtered.qzv ${metadata} asv_seqs_tax_filtered.qzv 09_filter_table_by_tax_output ${tax_tsv} asv_table_tax_filtered.biom asv_table_tax_filtered.tsv completecmd &> q2_filter_table_by_tax.log 2>&1
    """

}

process q2_filter_table_by_data {

    label 'qiime2_env'

    publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: '10_filter_table_by_data_output'
    publishDir "${params.outdir}/${params.filter_table_by_data_step}", mode: 'copy', pattern: '*.qz*'
    publishDir "${params.outdir}/${params.filter_table_by_data_step}", mode: 'copy', pattern: 'asv_table_filtered.biom'
    publishDir "${params.outdir}/${params.report_dirname}/99_completecmd", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_filtering_data -> "cmd/10_${task.process}_complete.sh" }

    input:
        path(asv_table)
        path(asv_seqs)
        path(metadata)
        path(tax_tsv)

    output:
        path('final_asv_table_filtered.qza'), emit: final_asv_table_filtered_qza
        path('final_asv_table_filtered.qzv')
        path('final_asv_seqs_filtered.qza'), emit: final_asv_seqs_filtered_qza
        path('final_asv_seqs_filtered.qzv')
        path('10_filter_table_by_data_output'), emit: filter_table_by_data_output
        path('10_filter_table_by_data_output/sequences.fasta'), emit: filter_table_by_data_seqs_fasta
        path('asv_table_filtered.biom'), emit: final_asv_table_filtered_biom
        path('asv_table_filtered.tsv'), emit: final_asv_table_filtered_tsv
        path('completecmd')

    script:
    """
    10_q2_filter_table_by_data.sh ${params.filter_by_id} ${params.list_sample_to_remove} ${asv_table} asv_table_samples_filtered.qza ${asv_seqs} asv_seqs_samples_filtered.qza ${params.filter_by_frequency} ${params.min_frequency_sample} ${params.min_frequency_asv} ${params.contingency_asv} asv_table_frequency_filtered.qza asv_seqs_frequency_filtered.qza final_asv_table_filtered.qzv ${metadata} final_asv_seqs_filtered.qzv 10_filter_table_by_data_output ${tax_tsv} asv_table_filtered.biom asv_table_filtered.tsv completecmd &> q2_filter_table_by_data.log 2>&1
    """

}

process q2_asv_phylogeny {

    label 'qiime2_env'
    label 'multithreads'

    publishDir "${params.outdir}/${params.asv_phylogeny_results}", mode: 'copy', pattern: 'asv_phylogeny.nwk'
    publishDir "${params.outdir}/${params.steps_data_dirname}", mode: 'copy', pattern: 'asv_phylogeny'
    publishDir "${params.outdir}/${params.report_dirname}/98_version", mode: 'copy', pattern : 'v_mafft.txt'
    publishDir "${params.outdir}/${params.report_dirname}/98_version", mode: 'copy', pattern : 'v_fasttree.txt'
    publishDir "${params.outdir}/${params.report_dirname}/99_completecmd", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_phylo -> "cmd/12_${task.process}_complete.sh" }

    input:
        path(asv_seqs)

    output:
        path('asv_phylogeny')
        path('asv_phylogeny.nwk'), emit: asv_phylogeny_nwk
        path('v_mafft.txt')
        path('v_fasttree.txt')
        path('completecmd')

    script:
    """
    12_q2_phylogeny.sh ${task.cpus} ${asv_seqs} asv_phylogeny asv_phylogeny.nwk completecmd >& q2_phylogeny.log 2>&1
    mafft --version > v_mafft.txt
    echo '2.1.11' > v_fasttree.txt
    """

}

process q2_ancombc {

    tag "${ancombc_formula}"
    label 'qiime2_env'

    publishDir "${params.outdir}/${params.ancombc_results}", mode: 'copy', pattern: '*_level_ancombc_*'
    publishDir "${params.outdir}/${params.ancombc_results}", mode: 'copy', pattern: 'heatmap_ancombc_*'
    publishDir "${params.outdir}/${params.ancombc_step}", mode: 'copy', pattern: '*.qz*'
    publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'completecmd', saveAs : { completecmd_ancombc -> "cmd/13_${task.process}_complete.sh" }

    input:
        path(asv_table)
        path(metadata)
        path(taxonomy)
        each(ancombc_formula)

    output:
        path('ancombc_output_*.qza')
        path('*_level_ancombc_*')
        path('ancombc_family_output_*.qza')
        path('ancombc_genus_output_*.qza')
        path('heatmap_ancombc_*')
        path('completecmd')

    script:
    def ref_level = params.use_custom_reference ? params.reference_level : "none"
    """
    13a_q2_ancombc.sh ${asv_table} ${metadata} ${params.use_custom_reference} ${ref_level} ${params.p_adj_method} ${params.max_iter} ${params.alpha} ${ancombc_formula} ancombc_output_${ancombc_formula}.qza asv_level_ancombc_${ancombc_formula} ${taxonomy} ancombc_table_family_${ancombc_formula}.qza ancombc_family_output_${ancombc_formula}.qza family_level_ancombc_${ancombc_formula} ancombc_table_genus_${ancombc_formula}.qza ancombc_genus_output_${ancombc_formula}.qza genus_level_ancombc_${ancombc_formula} completecmd &> q2_ancom.log 2>&1
    Rscript --vanilla ${baseDir}/bin/13b_ancombc_summary.R ${ancombc_formula} asv_level_ancombc_${ancombc_formula}/lfc_slice.csv asv_level_ancombc_${ancombc_formula}/q_val_slice.csv asv_level_ancombc_${ancombc_formula}_summary.tsv ${metadata} heatmap_ancombc_asv_level_${ancombc_formula}
    Rscript --vanilla ${baseDir}/bin/13b_ancombc_summary.R ${ancombc_formula} family_level_ancombc_${ancombc_formula}/lfc_slice.csv family_level_ancombc_${ancombc_formula}/q_val_slice.csv family_level_ancombc_${ancombc_formula}_summary.tsv ${metadata} heatmap_ancombc_family_level_${ancombc_formula}
    Rscript --vanilla ${baseDir}/bin/13b_ancombc_summary.R ${ancombc_formula} genus_level_ancombc_${ancombc_formula}/lfc_slice.csv genus_level_ancombc_${ancombc_formula}/q_val_slice.csv genus_level_ancombc_${ancombc_formula}_summary.tsv ${metadata} heatmap_ancombc_genus_level_${ancombc_formula}
    """

}
