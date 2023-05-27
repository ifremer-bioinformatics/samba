process picrust2 {

    label 'picrust2_env'
    label 'medRAM'

    publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: '14_PICRUSt2_predictions_output'
    publishDir "${params.outdir}/${params.report_dirname}/98_version", mode: 'copy', pattern: 'v_*.txt'
    publishDir "${params.outdir}/${params.report_dirname}/99_completecmd", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_picrust2 -> "cmd/14a_${task.process}_complete.sh" }

    input:
        path(asv_table_biom)
        path(asv_seqs_fasta)
        path(asv_table)
        path(metadata)
        each(var)

    output:
        path('14_PICRUSt2_predictions_output'), emit: picrust2_outdir
        path('14_PICRUSt2_predictions_output/RDA_plots/*predictions_RDA_*')
        path('v_*.txt')

    script:
    """
    14a_picrust2_process.sh ${asv_seqs_fasta} ${asv_table_biom} 14_PICRUSt2_predictions_output ${task.cpus} ${params.traits_db} ${params.nsti} ${params.hsp_method} ${params.min_reads} ${params.min_samples} ${params.top_level_mapfile} ${params.secondary_level_mapfile} completecmd >& picrust2_process.log 2>&1
    Rscript --vanilla ${baseDir}/bin/14b_picrust2_plots.R ${params.traits_db} ${asv_table} ${metadata} _predictions_RDA_ ${var} >& picrust2_plots.log 2>&1
    mkdir -p 14_PICRUSt2_predictions_output/RDA_plots
    mv *predictions_RDA_* 14_PICRUSt2_predictions_output/RDA_plots

    ## get statistics libraries version for report
    Rscript -e "library(svglite); write(x=as.character(packageVersion('svglite')), file='v_svglite.txt')"
    Rscript -e "library(ggord); write(x=as.character(packageVersion('ggord')), file='v_ggord.txt')"
    Rscript -e "library(compositions); write(x=as.character(packageVersion('compositions')), file='v_compositions.txt')"
    """

}
