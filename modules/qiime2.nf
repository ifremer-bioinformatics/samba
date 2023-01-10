process q2_import_data {

    label 'qiime2_env'

    publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: '03_import_output'
    publishDir "${params.outdir}/${params.import_step}", mode: 'copy', pattern: 'data.qz*'
    publishDir "${params.outdir}/${params.report_dirname}/98_version", mode: 'copy', pattern: 'v_*.txt'
    publishDir "${params.outdir}/${params.report_dirname}/99_completecmd", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_import -> "${task.process}_complete.sh" }

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
    publishDir "${params.outdir}/${params.report_dirname}/99_completecmd", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_cutadapt -> "cmd/${task.process}_complete.sh" }

    input:
        path(imported_data)

    output:
        path('trimmed_data.qza'), emit: trimmed_data
        path('trimmed_data.qzv')
        path('04_cutadapt_output')
        path('completecmd')

    script:
    """
    04_q2_cutadapt.sh ${params.singleEnd} ${task.cpus} ${imported_data} ${params.primerF} ${params.primerR} ${params.errorRate} trimmed_data.qza trimmed_data.qzv 04_cutadapt_output completecmd &> q2_cutadapt.log 2>&1
    """

}

process q2_dada2 {

    label 'qiime2_env'
    label 'multithreads'

    publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: '06_DADA2_output'
    publishDir "${params.outdir}/${params.dada2_step}", mode: 'copy', pattern: 'DADA2_*'
    publishDir "${params.outdir}/${params.report_dirname}/99_completecmd", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_dada2 -> "cmd/${task.process}_complete.sh" }

    input:
        path(dada2_input_data)
        path(metadata)

    output:
        path('DADA2_rep_seqs.qza'), emit: dada2_rep_seqs
        path('DADA2_rep_seqs.qzv')
        path('DADA2_table.qza'), emit: dada2_table
        path('DADA2_table.qzv')
        path('DADA2_process_stats.qz*')
        path('06_DADA2_output'), emit: dada2_outdir
        path('completecmd')

    script:
    """
    06_q2_dada2.sh ${params.singleEnd} ${dada2_input_data} ${params.FtrimLeft} ${params.RtrimLeft} ${params.FtruncLen} ${params.RtruncLen} ${params.truncQ} ${params.FmaxEE} ${params.RmaxEE} ${params.pooling_method} ${params.chimeras_method} ${task.cpus} DADA2_rep_seqs.qza DADA2_table.qza DADA2_process_stats.qza DADA2_process_stats.qzv DADA2_table.qzv ${metadata} DADA2_rep_seqs.qzv 06_DADA2_output completecmd &> q2_dada2.log 2>&1
"""

}

process q2_dbOTU3 {

    label 'qiime2_env'

    publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: '07_dbOTU3_output'
    publishDir "${params.outdir}/${params.dbotu3_step}", mode: 'copy', pattern: 'dbOTU3*'
    publishDir "${params.outdir}/${params.report_dirname}/99_completecmd", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_dbotu3 -> "cmd/${task.process}_complete.sh" }

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
        path('dbOTU3_output'), emit: dbotu3_outdir
        path('completecmd')

    script:
    """
    07_q2_dbotu3.sh ${asv_table} ${asv_seqs} ${params.gen_crit} ${params.abund_crit} ${params.pval_crit} dbOTU3_seqs.qza dbOTU3_table.qza dbOTU3_details.txt dbOTU3_table.qzv ${metadata} dbOTU3_seqs.qzv dbOTU3_output completecmd &> q2_dbotu3.log 2>&1
    """

}

process q2_assign_taxo {

    label 'qiime2_env'
    label 'highRAM'

    publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: '08_taxonomic_assignation_output'
    publishDir "${params.outdir}/${params.assign_taxo_step}", mode: 'copy', pattern: '*.qz*'
    publishDir "${params.outdir}/${params.assign_taxo_step}", mode: 'copy', pattern: 'ASV_tax_table.biom'
    publishDir "${params.outdir}/${params.report_dirname}/99_completecmd", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_assign_taxo -> "cmd/${task.process}_complete.sh" }

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
    publishDir "${params.outdir}/${params.report_dirname}/99_completecmd", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_filtering_tax -> "cmd/${task.process}_complete.sh" }

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
        path('asv_table_tax_filtered.biom'), emit: asv_table_tax_filtered_biom
        path('asv_table_tax_filtered.tsv'), emit: asv_table_tax_filtered_tsv
        path('completecmd')

    script:
    """
    09_q2_filter_table_by_tax.sh ${asv_table} ${tax_qza} ${params.tax_to_exclude} asv_table_tax_filtered.qza ${asv_seqs} asv_seqs_tax_filtered.qza asv_table_tax_filtered.qzv ${metadata} asv_seqs_tax_filtered.qzv 09_filter_table_by_tax_output ${tax_tsv} asv_table_tax_filtered.biom asv_table_tax_filtered.tsv completecmd &> q2_filter_table_by_tax.log 2>&1
    """

}

process q2_filter_table_by_data {

    label 'qiime2_env'

    publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: '10_filter_table_by_data_output'
    publishDir "${params.outdir}/${params.filter_table_by_data_step}", mode: 'copy', pattern: '*.qz*'
    publishDir "${params.outdir}/${params.filter_table_by_data_step}", mode: 'copy', pattern: 'biom'
    publishDir "${params.outdir}/${params.report_dirname}/99_completecmd", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_filtering_data -> "cmd/${task.process}_complete.sh" }

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
        path('asv_table_filtered.biom'), emit: final_asv_table_filtered_biom
        path('asv_table_filtered.tsv'), emit: final_asv_table_filtered_tsv
        path('completecmd')

    script:
    """
    10_q2_filter_table_by_data.sh ${params.filter_by_id} ${params.list_sample_to_remove} ${asv_table} asv_table_samples_filtered.qza ${asv_seqs} asv_seqs_samples_filtered.qza ${params.filter_by_frequency} ${params.min_frequency_sample} ${params.min_frequency_asv} ${params.contingency_asv} asv_table_frequency_filtered.qza asv_seqs_frequency_filtered.qza final_asv_table_filtered.qzv ${metadata} final_asv_seqs_filtered.qzv 10_filter_table_by_data_output ${tax_tsv} asv_table_filtered.biom asv_table_filtered.tsv completecmd &> q2_filter_table_by_data.log 2>&1
    """

}
