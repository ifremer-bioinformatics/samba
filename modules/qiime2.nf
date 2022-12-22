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
        path('06_DADA2_output')
        path('completecmd')

    script:
    """
    06_q2_dada2.sh ${params.singleEnd} ${dada2_input_data} ${params.FtrimLeft} ${params.RtrimLeft} ${params.FtruncLen} ${params.RtruncLen} ${params.truncQ} ${params.FmaxEE} ${params.RmaxEE} ${params.pooling_method} ${params.chimeras_method} ${task.cpus} DADA2_rep_seqs.qza DADA2_table.qza DADA2_process_stats.qza DADA2_process_stats.qzv DADA2_table.qzv ${metadata} DADA2_rep_seqs.qzv 06_DADA2_output completecmd &> q2_dada2.log 2>&1
"""

}
