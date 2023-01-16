process data_integrity {

    label 'biopython'

    publishDir "${params.outdir}/${params.data_integrity_results}", mode: 'copy', pattern: 'data_integrity.txt'
    publishDir "${params.outdir}/${params.data_integrity_results}", mode: 'copy', pattern: '*.sort'

    input:
        path(manifest)
        path(metadata)

    output:
        path('data_integrity.txt')
        path("q2_metadata_validated.tsv"), emit: final_metadata
        path("q2_manifest_validated.tsv"), emit: final_manifest

    script:
    def datatype = params.singleEnd ? "single" : "paired"
    def control = params.list_control_samples ? params.list_control_samples : "none"
    """
    02_data_integrity.py -e ${metadata} -a ${manifest} -p ${params.primer_filter} -s ${datatype} -t ${task.cpus} -c ${control} &> data_integrity.log 2>&1
    """
}
