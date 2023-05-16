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
    publishDir "${params.outdir}/${params.nanopore_r_results}/01_data", mode: 'copy', pattern: '*.rds'
    publishDir "${params.outdir}/${params.report_dirname}/98_version", mode: 'copy', pattern: 'v_*.txt'
    publishDir "${params.outdir}/${params.report_dirname}/99_completecmd", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_phyloseq -> "06_${task.process}_complete.sh" }

    input:
        path(asv_table_tsv)
        path(metadata)

    output:
        path('*.tsv')
        path('*.rds'), emit: phy_obj
        path('v_*.txt')
        path('completecmd')

    script:
    """
    Rscript --vanilla ${baseDir}/bin/NANOPORE_04_phyloseq.R phyloseq_ ${asv_table_tsv} ${metadata} count_table_for_stats_all_assignation.tsv count_table_for_stats_only_assigned.tsv ${params.db_name} &> stats_prepare_data.log 2&>1
    cp ${baseDir}/bin/NANOPORE_04_phyloseq.R completecmd

    ## get statistics libraries version for report
    Rscript -e "write(x=as.character(paste0(R.Version()[c('major','minor')], collapse = '.')), file='v_R.txt')"
    Rscript -e "library(dplyr); write(x=as.character(packageVersion('dplyr')), file='v_dplyr.txt')"
    Rscript -e "library(stringr); write(x=as.character(packageVersion('stringr')), file='v_stringr.txt')"
    Rscript -e "library(phyloseq); x=as.character(packageVersion('phyloseq')); write(x, file='v_phyloseq.txt')"
    """

}

process nanopore_alpha_diversity {

    tag "${var}"
    label 'R_env'

    publishDir "${params.outdir}/${params.nanopore_r_results}/02_analysis/01_alpha_diversity", mode: 'copy', pattern: 'alpha_div_index_values_*.txt'
    publishDir "${params.outdir}/${params.nanopore_r_results}/02_analysis/01_alpha_diversity/figures_all_assignation", mode: 'copy', pattern: 'alpha_div_bxp_all_assignation_*'
    publishDir "${params.outdir}/${params.nanopore_r_results}/02_analysis/01_alpha_diversity/figures_only_assigned", mode: 'copy', pattern: 'alpha_div_bxp_only_assigned_*'
    publishDir "${params.outdir}/${params.nanopore_r_results}/02_analysis/01_alpha_diversity/figures_all_assignation", mode: 'copy', pattern: 'rarefaction_curve_all_assignation_*'
    publishDir "${params.outdir}/${params.nanopore_r_results}/02_analysis/01_alpha_diversity/figures_only_assigned", mode: 'copy', pattern: 'rarefaction_curve_only_assigned_*'
    publishDir "${params.outdir}/${params.nanopore_r_results}/02_analysis/01_alpha_diversity/figures_all_assignation", mode: 'copy', pattern: 'abundance_table_all_assignation_*.tsv'
    publishDir "${params.outdir}/${params.nanopore_r_results}/02_analysis/01_alpha_diversity/figures_only_assigned", mode: 'copy', pattern: 'abundance_table_only_assigned_*.tsv'
    publishDir "${params.outdir}/${params.nanopore_r_results}/02_analysis/01_alpha_diversity/figures_all_assignation", mode: 'copy', pattern: 'pie_all_assignation_*'
    publishDir "${params.outdir}/${params.nanopore_r_results}/02_analysis/01_alpha_diversity/figures_only_assigned", mode: 'copy', pattern: 'pie_only_assigned_*'
    publishDir "${params.outdir}/${params.nanopore_r_results}/02_analysis/01_alpha_diversity/figures_all_assignation", mode: 'copy', pattern: 'sample_abundance_all_assignation_*.tsv'
    publishDir "${params.outdir}/${params.nanopore_r_results}/02_analysis/01_alpha_diversity/figures_only_assigned", mode: 'copy', pattern: 'sample_abundance_only_assigned_*.tsv'
    publishDir "${params.outdir}/${params.nanopore_r_results}/02_analysis/01_alpha_diversity/figures_all_assignation", mode: 'copy', pattern: 'barplot_all_assignation_*'
    publishDir "${params.outdir}/${params.nanopore_r_results}/02_analysis/01_alpha_diversity/figures_only_assigned", mode: 'copy', pattern: 'barplot_only_assigned_*'
    publishDir "${params.outdir}/${params.report_dirname}/98_version", mode: 'copy', pattern: 'v_*.txt'
    publishDir "${params.outdir}/${params.report_dirname}/99_completecmd", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_alpha_diversity -> "07_${task.process}_complete.sh" }

    input:
        path(phyloseq)
        each(var)

    output:
        path('alpha_div_index_values_*.txt')
        path('alpha_div_bxp_all_assignation_*')
        path('alpha_div_bxp_only_assigned_*')
        path('rarefaction_curve_all_assignation_*')
        path('rarefaction_curve_only_assigned_*')
        path('abundance_table_all_assignation_*.tsv')
        path('abundance_table_only_assigned_*.tsv')
        path('pie_all_assignation_*')
        path('pie_only_assigned_*')
        path('sample_abundance_all_assignation_*.tsv')
        path('sample_abundance_only_assigned_*.tsv')
        path('barplot_all_assignation_*')
        path('barplot_only_assigned_*')
        path('completecmd')
        path('v_*.txt')
        val('report_ok'), emit: report_ok

    script :
    """
    Rscript --vanilla ${baseDir}/bin/NANOPORE_05_alpha_diversity.R phyloseq_all_assignation alpha_div_index_values_all_assignation.txt ${var} alpha_div_bxp_all_assignation rarefaction_curve_all_assignation ${params.taxa_nb} abundance_table_all_assignation pie_all_assignation sample_abundance_all_assignation barplot_all_assignation ${params.db_name} > alpha_diversity_all_assignation.log 2&>1
    Rscript --vanilla ${baseDir}/bin/NANOPORE_05_alpha_diversity.R phyloseq_only_assigned alpha_div_index_values_only_assigned.txt ${var} alpha_div_bxp_only_assigned rarefaction_curve_only_assigned ${params.taxa_nb} abundance_table_only_assigned pie_only_assigned sample_abundance_only_assigned barplot_only_assigned ${params.db_name} > alpha_diversity_only_assigned.log 2&>1
    cp ${baseDir}/bin/NANOPORE_05_alpha_diversity.R completecmd
    touch report_ok

    ## get statistics libraries version for report
    Rscript -e "library(tidyr); write(x=as.character(packageVersion('tidyr')), file='v_tidyr.txt')"
    Rscript -e "library(ggplot2); write(x=as.character(packageVersion('ggplot2')), file='v_ggplot2.txt')"
    Rscript -e "library(rstatix); write(x=as.character(packageVersion('rstatix')), file='v_rstatix.txt')"
    Rscript -e "library(ggpubr); write(x=as.character(packageVersion('ggpubr')), file='v_ggpubr.txt')"
    Rscript -e "library(grid); write(x=as.character(packageVersion('grid')), file='v_grid.txt')"
    Rscript -e "library(vegan); x=as.character(packageVersion('vegan')); write(x, file='v_vegan.txt')"
    Rscript -e "library(RColorBrewer); x=as.character(packageVersion('RColorBrewer')); write(x, file='v_RColorBrewer.txt')"
    """

}
