process filter_contaminants {
    
    label 'qiime2_env'

    publishDir "${params.outdir}/${params.filter_contaminants_results}", mode: 'copy', pattern: 'decontaminated_ASV_table.tsv'
    publishDir "${params.outdir}/${params.filter_contaminants_results}", mode: 'copy', pattern: 'abundance_removed.txt'
    publishDir "${params.outdir}/${params.filter_contaminants_results}", mode: 'copy', pattern: 'ASV_removed.txt'
    publishDir "${params.outdir}/${params.filter_contaminants_results}", mode: 'copy', pattern: 'decontaminated_ASV.fasta'
    publishDir "${params.outdir}/${params.filter_contaminants_results}", mode: 'copy', pattern: 'filter_contaminants_export'
    publishDir "${params.outdir}/${params.filter_contaminants_step}", mode: 'copy', pattern: '*.qz*'
    publishDir "${params.outdir}/${params.filter_contaminants_step}", mode: 'copy', pattern: 'decontaminated_ASV_table.biom'
    publishDir "${params.outdir}/${params.report_dirname}/98_version", mode: 'copy', pattern: 'v_microdecon.txt'
    publishDir "${params.outdir}/${params.report_dirname}/99_completecmd", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_filter_contaminants -> "cmd/11_${task.process}_complete.sh" }
 
    input :
        path(asv_table)
        path(asv_seqs)
        path(metadata)

    output :
        path('decontaminated_ASV_table.tsv'), emit: decontam_ASV_table_tsv
        path('abundance_removed.txt')
        path('ASV_removed.txt')
        path('decontaminated_ASV_table.biom'), emit: decontam_ASV_table_biom
        path('decontaminated_ASV_table.qza'), emit: decontam_ASV_table_qza
        path('decontaminated_ASV_table.qzv')
        path('decontaminated_ASV.fasta'), emit: decontam_ASV_seqs_fasta
        path('decontaminated_ASV_seqs.qza'), emit: decontam_ASV_seqs_qza
        path('decontaminated_ASV_seqs.qzv')
        path('filter_contaminants_export')
        path('completecmd')
        path('v_microdecon.txt')

    script :
    """
    sed '1d' ${asv_table} > microdecon_table
    sed -i 's/#OTU ID/ASV_ID/g' microdecon_table
    sed -i "s/'//g" microdecon_table
    Rscript --vanilla ${baseDir}/bin/11a_microdecon.R microdecon_table ${params.list_control_samples} decontaminated_ASV_table.tsv abundance_removed.txt ASV_removed.txt &> microdecon.log 2>&1
    cp ${baseDir}/bin/11a_microdecon.R completecmd
    Rscript -e "write(x=as.character(packageVersion('microDecon')), file='v_microdecon.txt')"
    11b_microdecon_output.sh decontaminated_ASV_table.tsv decontaminated_ASV_table.biom decontaminated_ASV_table.qza ${asv_seqs} decontaminated_ASV.fasta decontaminated_ASV_seqs.qza decontaminated_ASV_table.qzv ${metadata} filter_contaminants_export decontaminated_ASV_seqs.qzv completecmd &> microdecon-to_qiime2.log 2>&1
    """
}
