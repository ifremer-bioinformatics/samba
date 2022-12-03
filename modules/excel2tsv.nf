process excel2tsv {

    label 'biopython'

    publishDir "${params.outdir}/${params.excel2tsv_results}", mode: 'copy', pattern: 'q2_manifest.sort.tsv'
    publishDir "${params.outdir}/${params.excel2tsv_results}", mode: 'copy', pattern: 'q2_metadata.sort.tsv'
    publishDir "${params.outdir}/${params.excel2tsv_results}", mode: 'copy', pattern: 'excel2tsv.log'

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
    01_excel2tsv.py -x ${xls} -s ${datatype} &> excel2tsv.log 2>&1
    """

}

process addpath_testdata {

    publishDir "${params.outdir}/${params.excel2tsv_dirname}", mode: 'copy', pattern: 'manifest.tsv'

    input:
        path(manifest)

    output:
        path("manifest.tsv"), emit: manifest

    script:
    """
    01_add_path_test_data.sh ${manifest} ${baseDir} manifest.tsv completecmd &> add_path_test_data.log 2>&1
    """

}
