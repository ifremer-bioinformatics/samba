process format_final_outputs {
    
    label 'qiime2_env'

    publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'final_asv_table.tsv'
    publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'final_asv_sequences.fasta'
 
    input :
        path(asv_table_tsv)
        path(asv_seqs_fasta)

    output :
        path('final_asv_table.tsv'), emit: final_asv_table
        path('final_asv_sequences.fasta')

    script :
    """
    if "${params.filter_contaminants_enable}"
    then
        cp ${asv_table_tsv} final_asv_table.tsv
        sed -i "s/'//g" final_asv_table.tsv
        sed -i 's/; /;/g' final_asv_table.tsv
    else
        sed '1d' ${asv_table_tsv} > final_asv_table.tsv
        sed -i 's/#OTU ID/ASV_ID/g' final_asv_table.tsv
        sed -i "s/'//g" final_asv_table.tsv
        sed -i 's/; /;/g' final_asv_table.tsv
    fi
    cp ${asv_seqs_fasta} final_asv_sequences.fasta
    """

}
