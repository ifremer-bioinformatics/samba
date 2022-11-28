process excel2tsv {

    label 'biopython_env'

    publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: 'q2_manifest.sort.tsv'
    publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: 'q2_metadata.sort.tsv'
    publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: 'excel2tsv.log'

    input:
        path(xls)
        val(ready)

    output:
        path("q2_metadata.sort.tsv"), emit: metadata_xls
        path("q2_manifest.sort.tsv"), emit: manifest_xls
        path("excel2tsv.log")

    script:
    def datatype = params.singleEnd ? "single" : "paired"
    """
    excel2tsv.py -x ${xls} -s ${datatype} &> excel2tsv.log 2>&1
    """

}
