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
    Rscript --vanilla ${baseDir}/bin/15_create_phyloseq_obj.R ${asv_table_tsv} asv_table_for_stats.tsv ${metadata} metadata_for_stats.tsv ${tree} phyloseq.rds ${params.db_name} &> stats_prepare_data.log 2&>1
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
    publishDir "${params.outdir}/${params.report_dirname}/99_completecmd", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_phyloseq -> "06a_${task.process}_complete.sh" }

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
    Rscript --vanilla ${baseDir}/bin/NANOPORE_04a_phyloseq.R phyloseq_ ${asv_table_tsv} ${metadata} count_table_for_stats_all_assignation.tsv count_table_for_stats_only_assigned.tsv ${params.db_name} ${tax_to_filter} &> stats_prepare_data.log 2&>1
    cp ${baseDir}/bin/NANOPORE_04a_phyloseq.R completecmd

    ## get statistics libraries version for report
    Rscript -e "write(x=as.character(paste0(R.Version()[c('major','minor')], collapse = '.')), file='v_R.txt')"
    Rscript -e "library(dplyr); write(x=as.character(packageVersion('dplyr')), file='v_dplyr.txt')"
    Rscript -e "library(stringr); write(x=as.character(packageVersion('stringr')), file='v_stringr.txt')"
    Rscript -e "library(phyloseq); x=as.character(packageVersion('phyloseq')); write(x, file='v_phyloseq.txt')"
    """

}

process agglomerate_phyloseq {

    tag "${tax_level}"
    label 'R_env'

    publishDir "${params.outdir}/${params.nanopore_r_results}/01_data", mode: 'copy', pattern: '*.tsv'
    publishDir "${params.outdir}/${params.nanopore_r_results}/01_data", mode: 'copy', pattern: '*.rds'
    publishDir "${params.outdir}/${params.report_dirname}/99_completecmd", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_agglomerate_phyloseq -> "06b_${task.process}_complete.sh" }

    input:
        path(phyloseq)
        each(tax_level)

    output:
        path('*.tsv')
        path('*.rds'), emit: phy_obj_taxlevel
        path('completecmd')

    script:
    """
    Rscript --vanilla ${baseDir}/bin/NANOPORE_04b_agglomerate_taxlevel.R phyloseq_all_assignation.rds phyloseq_only_assigned.rds ${tax_level} &> agglomerate_phyloseq_taxlevel.log 2&>1
    cp ${baseDir}/bin/NANOPORE_04b_agglomerate_taxlevel.R completecmd
    """

}

process illumina_alpha_diversity {

    tag "${var}"
    label 'R_env'

    publishDir "${params.outdir}/${params.r_results}/02_analysis/01_alpha_diversity", mode: 'copy', pattern: 'alpha_div_index_values.txt'
    publishDir "${params.outdir}/${params.r_results}/02_analysis/01_alpha_diversity/${var}", mode: 'copy', pattern: 'alpha_div_bxp_*'
    publishDir "${params.outdir}/${params.r_results}/02_analysis/01_alpha_diversity/${var}/rarefaction_curve", mode: 'copy', pattern: 'rarefaction_curve_*'
    publishDir "${params.outdir}/${params.r_results}/02_analysis/01_alpha_diversity/${var}/barplots", mode: 'copy', pattern: 'barplot_*'
    publishDir "${params.outdir}/${params.r_results}/02_analysis/01_alpha_diversity/${var}/pies", mode: 'copy', pattern: 'pie_*'
    publishDir "${params.outdir}/${params.report_dirname}/98_version", mode: 'copy', pattern: 'v_*.txt'
    publishDir "${params.outdir}/${params.report_dirname}/99_completecmd", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_alpha_diversity -> "16_${task.process}_complete.sh" }

    input:
        path(phyloseq)
        each(var)

    output:
        path('alpha_div_index_values.txt')
        path('alpha_div_bxp_*')
        path('rarefaction_curve_*')
        path('pie_*')
        path('barplot_*')
        path('completecmd')
        path('v_*.txt')
        val('report_ok'), emit: report_ok

    script :
    """
    Rscript --vanilla ${baseDir}/bin/16_alpha_diversity.R ${phyloseq} alpha_div_index_values.txt ${var} alpha_div_bxp rarefaction_curve ${params.taxa_nb} ${params.db_name} > alpha_diversity.log 2&>1
    cp ${baseDir}/bin/16_alpha_diversity.R completecmd
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
        path(phyloseq_taxlevel)
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

process illumina_beta_diversity {

    tag "${var}_${normalisation}"
    label 'R_env'

    publishDir "${params.outdir}/${params.r_results}/01_data", mode: 'copy', pattern: 'asv_table_for_stats_*'
    publishDir "${params.outdir}/${params.r_results}/02_analysis/02_beta_diversity/${normalisation}/${var}", mode: 'copy', pattern: '*_pairwiseAdonis_result_*'
    publishDir "${params.outdir}/${params.r_results}/02_analysis/02_beta_diversity/${normalisation}/${var}", mode: 'copy', pattern: '*_permanova_result_*'
    publishDir "${params.outdir}/${params.r_results}/02_analysis/02_beta_diversity/${normalisation}/${var}/NMDS", mode: 'copy', pattern: 'ordination_NMDS_*'
    publishDir "${params.outdir}/${params.r_results}/02_analysis/02_beta_diversity/${normalisation}/${var}/PCoA", mode: 'copy', pattern: 'ordination_PCoA_*'
    publishDir "${params.outdir}/${params.r_results}/02_analysis/02_beta_diversity/${normalisation}", mode: 'copy', pattern: 'permanova_table_*.png'
    publishDir "${params.outdir}/${params.report_dirname}/99_completecmd", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_ordination_plots -> "17_${task.process}_complete.sh" }

    input:
        path(phyloseq)
        each(normalisation)
        each(var)

    output:
        path('asv_table_for_stats_*')
        path('ordination_*'), emit: ordi_plots
        path('*_permanova_result_*')
        path('*_pairwiseAdonis_result_*')
        path('permanova_table_*.png')
        path('completecmd')
        val('report_ok'), emit: report_ok

    script:
    """
    Rscript --vanilla ${baseDir}/bin/17_beta_diversity.R ${phyloseq} ${normalisation} ${var} &> beta_diversity.log 2&>1
    cp ${baseDir}/bin/17_beta_diversity.R completecmd
    touch report_ok
    """

}

process intersecting_sets {

    tag "${var}"
    label 'R_env'

    publishDir "${params.outdir}/${params.r_results}/02_analysis/03_intersecting_sets/${var}", mode: 'copy', pattern: '*_UpSetR_*.png'
    publishDir "${params.outdir}/${params.report_dirname}/99_completecmd", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_pie_class -> "18_${task.process}_complete.sh" }

    input:
        path(phyloseq)
        each(var)

    output :
        path('*_UpSetR_*.png')
        path('completecmd')
        val('report_ok'), emit: report_ok

    script:
    """
    Rscript --vanilla ${baseDir}/bin/18_intersecting_sets.R ${phyloseq} ${var} ${params.db_name} &> UpSetR.log 2&>1
    cp ${baseDir}/bin/18_intersecting_sets.R completecmd
    touch report_ok
    """

}
