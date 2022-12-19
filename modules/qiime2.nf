process import_data {

    label 'qiime2_env1cpu'

    publishDir "${params.outdir}/${params.import_results}", mode: 'copy', pattern: '*_output'
    publishDir "${params.outdir}/${params.import_step}", mode: 'copy', pattern: 'data.qz*'
    publishDir "${params.outdir}/${params.report_dirname}/98_version", mode: 'copy', pattern: 'v_*.txt'
    publishDir "${params.outdir}/${params.report_dirname}/99_completecmd", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_import -> "${task.process}_complete.sh" }

    input:
        path(final_manifest)

    output:
        path('data.qza'), emit: imported_data
        path('data.qzv')
        path('import_output/*')
        path('completecmd')
        path('v_qiime2.txt')

    script:
    """
    03_q2_import.sh ${params.singleEnd} ${final_manifest} data.qza data.qzv import_output completecmd &> q2_import.log 2>&1
    qiime --version|grep 'q2cli'|cut -d' ' -f3 > v_qiime2.txt
    """
}
