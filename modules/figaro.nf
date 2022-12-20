process figaro {

    label 'figaro_env'

    publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: '05_figaro_output'
    publishDir "${params.outdir}/${params.report_dirname}/99_completecmd", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_figaro -> "${task.process}_complete.sh" }

    input:
        val(ready)

    output:
        path('05_figaro_output')
        path('completecmd')

    script:
    """
    05_figaro.sh ${params.raw_data_dir} ${params.amplicon_length} ${params.primerF} ${params.primerR} 05_figaro_output completecmd &> figaro.log 2>&1
    """

}
