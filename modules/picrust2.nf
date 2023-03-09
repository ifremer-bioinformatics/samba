process picrust2 {

    label 'picrust2_env'
    label 'medRAM'

    publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: '14_PICRUSt2_predictions_output'
    publishDir "${params.outdir}/${params.report_dirname}/99_completecmd", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_picrust2 -> "cmd/14a_${task.process}_complete.sh" }

    input:
        path(asv_table_biom)
        path(asv_seqs_fasta)

    output:
        path('14_PICRUSt2_predictions_output'), emit: picrust2_outdir

    script:
    """
    14a_picrust2_process.sh ${asv_seqs_fasta} ${asv_table_biom} 14_PICRUSt2_predictions_output ${task.cpus} ${params.traits_db} ${params.nsti} ${params.hsp_method} ${params.min_reads} ${params.min_samples} completecmd >& picrust2_process.log 2>&1
    """

}
