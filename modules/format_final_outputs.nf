process format_final_outputs {
    
    label 'qiime2_env'

    publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'final_asv_table.tsv'
    publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'final_asv_sequences.fasta'
    publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'final_asv_table.biom'
    publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'final_asv_table.qza'
    publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'final_asv_sequences.qza'
 
    input :
        path(asv_table_tsv)
        path(asv_seqs_fasta)
        path(asv_table_biom)
        path(asv_table_qza)
        path(asv_seqs_qza)

    output :
        path('final_asv_table.tsv'), emit: final_asv_table
        path('final_asv_sequences.fasta')
        path('final_asv_table.biom')
        path('final_asv_table.qza')
        path('final_asv_sequences.qza')

    script :
    """
    cp ${asv_table_tsv} final_asv_table.tsv
    sed -i "s/'//g" final_asv_table.tsv
    sed -i 's/; /;/g' final_asv_table.tsv
    cp ${asv_seqs_fasta} final_asv_sequences.fasta
    cp ${asv_table_biom} final_asv_table.biom
    cp ${asv_table_qza} final_asv_table.qza
    cp ${asv_seqs_qza} final_asv_sequences.qza
    """

}
