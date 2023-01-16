process swarm_clustering_processing {

    label 'swarm_env'
    label 'multithreads'

    publishDir "${params.outdir}/${params.report_dirname}/98_version", mode: 'copy', pattern: 'v_swarm.txt'
    publishDir "${params.outdir}/${params.report_dirname}/99_completecmd", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_swarm_clustering_processing -> "cmd/07b_${task.process}_complete.sh" }

    input:
        path(asv_seqs_fasta)

    output:
        path('final_asv_swarm_seqs_clusters.fasta'), emit: asv_swarm_seqs_clusters_fasta
        path('asv_swarm_cluster_list.tsv'), emit: asv_swarm_cluster_list_tsv

    script:
    """
    07b_swarm_clustering_processing.sh ${task.cpus} ${asv_seqs_fasta} asv_swarm_seqs_clusters.fasta tmp_swarm_cluster_list.txt final_asv_swarm_seqs_clusters.fasta completecmd >& swarm_clustering_processing.log 2>&1
    cat tmp_swarm_cluster_list.txt | awk '{ print length, \$0 }' | sort -nr -s | cut -d ' ' -f2- > asv_swarm_cluster_list.tsv
    echo '3.1.3' > v_swarm.txt
    """

}

process swarm_clustering_format_output {

    label 'qiime2_env'

    publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: '07_swarm_clustering'
    publishDir "${params.outdir}/${params.swarm_clustering_step}", mode: 'copy', pattern: '*.qz*'
    publishDir "${params.outdir}/${params.report_dirname}/99_completecmd", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_swarm_clustering_format_output -> "cmd/07cd_${task.process}_complete.sh" }

    input:
        path(asv_table)
        path(swarm_cluster_list)
        path(swarm_asv_seqs_fa)
        path(metadata)

    output:
        path('swarm_cluster_identification.tsv')
        path('asv_table_swarm_clustered.tsv'), emit: asv_table_swarm_clustered_tsv
        path('swarm_asv_table.qza'), emit: swarm_asv_table
        path('swarm_asv_table.qzv')
        path('07_swarm_clustering'), emit: swarm_asv_outdir
        path('07_swarm_clustering/sequences.fasta'), emit: swarm_asv_seqs_fasta
        path('swarm_asv_seqs.qza'), emit: swarm_asv_seqs
        path('swarm_asv_seqs.qzv')
        path('completecmd')

    script:
    """
    Rscript --vanilla ${baseDir}/bin/07c_swarm_clustering_table.R ${asv_table} ${swarm_cluster_list} swarm_cluster_identification.tsv asv_table_swarm_clustered.tsv >& swarm_clustering_table.log 2>&1
    cp ${baseDir}/bin/07c_swarm_clustering_table.R completecmd
    07d_swarm_clustering_output.sh asv_table_swarm_clustered.tsv asv_table_swarm_clustered.biom swarm_asv_table.qza ${swarm_asv_seqs_fa} swarm_asv_seqs.qza swarm_asv_table.qzv ${metadata} 07_swarm_clustering swarm_asv_seqs.qzv completecmd >& swarm_clustering_q2format.log 2>&1
    cp swarm_cluster_identification.tsv 07_swarm_clustering/
    cp asv_table_swarm_clustered.tsv 07_swarm_clustering/
    """

}
