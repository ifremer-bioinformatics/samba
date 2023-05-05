process figaro {

    label 'figaro_env'

    publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: '05_figaro_output'
    publishDir "${params.outdir}/${params.report_dirname}/98_version", mode: 'copy', pattern: 'v_figaro.txt'
    publishDir "${params.outdir}/${params.report_dirname}/99_completecmd", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_figaro -> "05_${task.process}_complete.sh" }

    input:
        val(ready)

    output:
        path('05_figaro_output')
        path('figaro.csv'), emit: figaro_csv
        path('v_figaro.txt')
        path('completecmd')

    script:
    """
    05a_figaro.sh ${params.raw_data_dir} ${params.amplicon_length} ${params.primerF} ${params.primerR} 05_figaro_output completecmd &> figaro.log 2>&1
    echo "1.1.2" > v_figaro.txt
    05b_figaro_json2csv.py -f 05_figaro_output/trimParameters.json &> parse_figaro_json.log 2>&1
    cp figaro.csv 05_figaro_output/
    """

}
