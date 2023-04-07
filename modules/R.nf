process create_phyloseq {

    label 'R_env'

    publishDir "${params.outdir}/${params.r_results}/01_data", mode: 'copy', pattern: '*.tsv'
    publishDir "${params.outdir}/${params.r_results}/01_data", mode: 'copy', pattern: 'phyloseq.rds'
    publishDir "${params.outdir}/${params.report_dirname}/98_version", mode: 'copy', pattern: 'v_*.txt'
    publishDir "${params.outdir}/${params.report_dirname}/99_completecmd", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_phyloseq -> "15_${task.process}_complete.sh" }

    input:
        path(asv_table_tsv)
        path(metadata)
        path(tree)

    output:
        path('asv_table_for_stats.tsv')
        path('metadata_for_stats.tsv')
        path('phyloseq.rds'), emit: phyloseq
        path('v_*.txt')
        path('completecmd')

    script:
    """
    Rscript --vanilla ${baseDir}/bin/15_create_phyloseq_obj.R ${asv_table_tsv} asv_table_for_stats.tsv ${metadata} metadata_for_stats.tsv ${tree} phyloseq.rds &> stats_prepare_data.log 2&>1
    cp ${baseDir}/bin/15_create_phyloseq_obj.R completecmd

    ## get statistics libraries version for report
    Rscript -e "write(x=as.character(paste0(R.Version()[c('major','minor')], collapse = '.')), file='v_R.txt')"
    Rscript -e "library(dplyr); write(x=as.character(packageVersion('dplyr')), file='v_dplyr.txt')"
    Rscript -e "library(stringr); write(x=as.character(packageVersion('stringr')), file='v_stringr.txt')"
    Rscript -e "library(phyloseq); x=as.character(packageVersion('phyloseq')); write(x, file='v_phyloseq.txt')"
    """

}

process nanopore_phyloseq_obj {

    label 'R_env'

    publishDir "${params.outdir}/${params.nanopore_r_results}/01_data", mode: 'copy', pattern: '*.tsv'
    publishDir "${params.outdir}/${params.nanopore_r_results}/01_data", mode: 'copy', pattern: 'phyloseq.rds'
    publishDir "${params.outdir}/${params.report_dirname}/99_completecmd", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_phyloseq -> "06_${task.process}_complete.sh" }

    input:
        path(asv_table_tsv)
        path(metadata)

    output:
        path('count_table_for_stats.tsv')
        path('metadata_for_stats.tsv')
        path('phyloseq.rds'), emit: phyloseq
        path('completecmd')

    script:
    """
    Rscript --vanilla ${baseDir}/bin/NANOPORE_04_phyloseq.R phyloseq.rds ${asv_table_tsv} ${metadata} ${params.tax_rank} ${params.kingdom} count_table_for_stats.tsv metadata_for_stats.tsv &> stats_prepare_data.log 2&>1
    cp ${baseDir}/bin/NANOPORE_04_phyloseq.R completecmd
    """

}
