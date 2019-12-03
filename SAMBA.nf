#!/usr/bin/env nextflow
println "Workflow for project : $params.projectName"
println "Workflow description : $workflow.manifest.description"
println "Workflow gitLab URL : $workflow.manifest.homePage"
println "Workflow authors : $workflow.manifest.author"
println "Workflow source code : $workflow.projectDir"
println "Cmd line: $workflow.commandLine"
println "Workflow working/temp directory : $workflow.workDir"
println "Workflow output/publish directory : $params.outdir"
println "Workflow configuration file : $workflow.configFiles"
println "Manifest file : $params.inmanifest"
println "Metadata file : $params.inmetadata"

Channel.fromPath(params.inmanifest, checkIfExists:true).into { manifest ; manifest4integrity }
Channel.fromPath(params.inmetadata, checkIfExists:true).into { metadata; metadata4stats ; metadata4integrity }

// IF NOT STATS ONLY, PERFORM QIIME STEPS
if(!params.stats_only){

    /* Check data integrity */
    
    process data_integrity {
        publishDir "${params.outdir}/${params.data_integrity_dirname}", mode: 'copy', pattern: 'data_integrity.csv'

    input :
        file manifest from manifest4integrity
        file metadata from metadata4integrity

    output :
        file 'verifications.ok' optional true into ckeck_ok
        file 'verifications.bad' optional true into check_bad
        file 'data_integrity.csv' into data_integrity_csv
    
    //Run only if process is activated in params.config file
    when :
        params.data_integrity_enable

    script :
    """
    ${baseDir}/lib/data_integrity.sh ${manifest} ${metadata} ${params.data_integrity.primerF} ${params.data_integrity.primerR} data_integrity.csv verifications.ok verifications.bad ${params.data_integrity.barcode} ${params.data_integrity.sampleid_column_name} ${params.data_integrity.R1_files_column_name} ${params.data_integrity.R2_files_column_name} ${params.data_integrity.barcode_filter} ${params.data_integrity.primer_filter} > data_integrity.log 2>&1
    if test -f "verifications.bad"; then
        echo "Data integrity process not satisfied, check ${params.outdir}/${params.data_integrity_dirname}/data_integrity.csv file"
        mkdir -p ${params.outdir}/${params.data_integrity_dirname}
        cp data_integrity.csv ${params.outdir}/${params.data_integrity_dirname}/.
        exit 1
     fi

    """
   }

    /* Import metabarcode data */
    
    process q2_import {
    
        beforeScript "${params.qiime_env}"
        publishDir "${params.outdir}/${params.import_dirname}", mode: 'copy', pattern: 'data.qz*'
        publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: '*_output'
        publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_import -> "cmd/${task.process}_complete.sh" }
    
        input : 
            file manifest from manifest
            file check_ok from ckeck_ok
    
        output : 
            file 'data.qza' into imported_data
            file 'data.qzv' into imported_visu
            file 'import_output' into imported_summary
            file 'completecmd' into complete_cmd_import

        //Run only if process is activated in params.config file
        when :
        params.qiime_import_enable
    
        script :
        """
        ${baseDir}/lib/q2_import.sh ${manifest} data.qza data.qzv import_output completecmd > q2_import.log 2>&1
        """
    }
    
    /* Trim metabarcode data with cutadapt */
    process q2_cutadapt {
    
        beforeScript "${params.qiime_env}"
        publishDir "${params.outdir}/${params.trimmed_dirname}", mode: 'copy', pattern: 'data*.qz*'
        publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: '*_output'
        publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_cutadapt -> "cmd/${task.process}_complete.sh" }
        
        input : 
            file imported_data from imported_data
    
        output :
            file 'data_trimmed.qza' into trimmed_data
            file 'data_trimmed.qzv' into trimmed_visu
            file 'trimmed_output' into trimmed_summary
            file 'completecmd' into complete_cmd_cutadapt
    
        //Run only if process is activated in params.config file
        when :
        params.cutadapt_enable
    
        script :
        """
        ${baseDir}/lib/q2_cutadapt.sh ${task.cpus} ${imported_data} ${params.cutadapt.primerF} ${params.cutadapt.primerR} ${params.cutadapt.errorRate} ${params.cutadapt.overlap} data_trimmed.qza data_trimmed.qzv trimmed_output completecmd > q2_cutadapt.log 2>&1
        """
    }

    /* Run dada2 */
    process q2_dada2 {
    
        beforeScript "${params.qiime_env}"
        publishDir "${params.outdir}/${params.dada2_dirname}", mode: 'copy', pattern: '*.qz*'
        publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: '*_output'
        publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_dada2 -> "cmd/${task.process}_complete.sh" }
    
        input : 
            file trimmed_data from trimmed_data
            file metadata from metadata 
    
        output :
            file 'rep_seqs.qza' into data_repseqs
            file 'rep_seqs.qzv' into visu_repseps
            file 'table.qza' into data_table
            file 'table.qzv' into visu_table
            file 'stats.qza' into stats_table
            file 'stats.qzv' into visu_stats
            file 'dada2_output' into dada2_summary
            file 'completecmd' into complete_cmd_dada2
    
        //Run only if process is activated in params.config file
        when :
        params.dada2_enable
    
        script :
        """
        ${baseDir}/lib/q2_dada2.sh ${trimmed_data} ${metadata} rep_seqs.qza rep_seqs.qzv table.qza table.qzv stats.qza stats.qzv dada2_output ${params.dada2.trim3F} ${params.dada2.trim3R} ${params.dada2.trunclenF} ${params.dada2.trunclenR} ${params.dada2.maxee_f} ${params.dada2.maxee_r} ${params.dada2.minqual} ${params.dada2.chimeras} ${task.cpus} completecmd > q2_dada2.log 2>&1
        """
    }

data_repseqs.into { repseqs_taxo ; repseqs_phylo }

    /* Run taxonomy assignment */

    process q2_taxonomy {
    
        beforeScript "${params.qiime_env}"
        publishDir "${params.outdir}/${params.taxo_dirname}", mode: 'copy', pattern: '*.qz*'
        publishDir "${params.outdir}/${params.taxo_dirname}", mode: 'copy', pattern: '*.tsv*'
        publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: '*_output'
        publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: 'Final_ASV_table*'
        publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_taxo -> "cmd/${task.process}_complete.sh" }
    
        input :
            file data_repseqs from repseqs_taxo
            file dada2_summary from dada2_summary
    
        output :
            file 'taxonomy.qza' into data_taxonomy
            file 'taxonomy.qzv' into visu_taxonomy
            file 'ASV_taxonomy.tsv' into taxonomy_tsv
            file 'taxo_output' into taxo_summary
            file 'Final_ASV_table_with_taxonomy.biom' into biom
            file 'Final_ASV_table_with_taxonomy.tsv' into biom_tsv
            file 'taxonomic_database.qza' optional true into trained_database
            file 'db_seqs_amplicons.qza' optional true into db_seqs_filtered
            file 'completecmd' into complete_cmd_taxo
    
        //Run only if process is activated in params.config file
        when :
        params.taxo_enable
    
        script :
        """
        ${baseDir}/lib/q2_taxo.sh ${task.cpus} ${params.taxo.db_seqs} ${params.taxo.db_tax} ${params.taxo.database} ${params.taxo.extract_db} ${params.cutadapt.primerF} ${params.cutadapt.primerR} ${params.taxo.confidence} ${data_repseqs} taxonomy.qza taxonomy.qzv taxo_output ASV_taxonomy.tsv ${dada2_summary} Final_ASV_table_with_taxonomy.biom Final_ASV_table_with_taxonomy.tsv taxonomic_database.qza db_seqs_amplicons.qza completecmd > q2_taxo.log 2>&1
        """ 
    }

    /* Run phylogeny construction */

    process q2_phylogeny {
        beforeScript "${params.qiime_env}"
        publishDir "${params.outdir}/${params.phylogeny_dirname}", mode: 'copy', pattern: '*.qza'
        publishDir "${params.outdir}/${params.phylogeny_dirname}", mode: 'copy', pattern: '*.txt'
        publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_phylo -> "cmd/${task.process}_complete.sh" }

        input :
            file repseqs from repseqs_phylo

        output :
            file 'aligned_repseq.qza' into aligned_repseq
            file 'masked-aligned_repseq.qza' into masked_aligned
            file 'tree.qza' into tree
            file 'model.txt' into model
            file 'tree_bestmodel.qza' into tree_bestmodel
            file 'tree_log' into tree_bestmodel_log
            file 'completecmd' into complete_cmd_phylogeny

        //Run only if process is activated in params.config file
        when :
        params.phylogeny_enable

        script :
        """
        ${baseDir}/lib/q2_phylogeny.sh repseqs aligned_repseq masked_aligned tree model ${params.phylogeny.alrt} ${params.phylogeny.bootstrap} tree_bestmodel tree_log ${task.cpus} completecmd > q2_phylogeny.log 2>&1
        """
    }
}

if(params.stats_only){
    
    //IF OTU TABLE ALREADY CREATED AND STAT ONLY STEPS NEEDED
    
    //Set biom_tsv path in params.conf
    Channel.fromPath(params.inasv_table, checkIfExists:true).set { biom_tsv }
    println "Input ASV table used for statistics steps : $params.inasv_table"
}

process prepare_data_for_stats {

    beforeScript "${params.r_stats_env}"

    publishDir "${params.outdir}/${params.report_dirname}/R/DATA", mode: 'copy', pattern : '*.tsv'
    publishDir "${params.outdir}/${params.report_dirname}/R/DATA", mode: 'copy', pattern : '*.rds'
    publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_prepare_stats -> "cmd/${task.process}_complete.sh" }
    
    input :
        file metadata from metadata4stats
        file biom_tsv from biom_tsv

    output :
        file 'ASV_table_with_taxo_for_stats.tsv' into biom_tsv_stats
        file 'metadata_stats.tsv' into metadata_stats
        file 'phyloseq.rds' into phyloseq_rds
        file 'completecmd' into complete_cmd_prepare_stats
 
    when :
    params.prepare_data_for_stats_enable

    script :
    """
    ${baseDir}/lib/prepare_data_for_stats.sh ${metadata} ${biom_tsv} ASV_table_with_taxo_for_stats.tsv metadata_stats.tsv completecmd > stats_prepare_data.log 2&>1
    Rscript --vanilla ${baseDir}/lib/create_phyloseq_obj.R phyloseq.rds ASV_table_with_taxo_for_stats.tsv metadata_stats.tsv >> stats_prepare_data.log 2&>1 
    """
}

//Duplicate channels needed in several processes
metadata_stats.into { metadata_beta ; metadata_beta_rarefied ; metadata_beta_deseq2 ; metadata_beta_css }
phyloseq_rds.into { phyloseq_rds_alpha ; phyloseq_rds_beta ; phyloseq_rds_beta_rarefied ; phyloseq_rds_beta_deseq2 ; phyloseq_rds_beta_css }

process stats_alpha {

    beforeScript "${params.r_stats_env}"

    publishDir "${params.outdir}/${params.report_dirname}/R/SCRIPT", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_alpha -> "${task.process}.R" }
    publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/alpha_diversity", mode: 'copy', pattern : '*.svg'
    
    input :
        file phyloseq_rds from phyloseq_rds_alpha

    output :
        file 'alpha_div_plots_*.svg' into alpha_div_plots
        file 'barplot_relabund_phylum_*.svg' into barplot_relabund_phylum
        file 'barplot_relabund_family_*.svg' into barplot_relabund_family
        file 'barplot_relabund_genus_*.svg' into barplot_relabund_genus
        file 'completecmd' into complete_cmd_alpha

    //Run only if process is activated in params.config file
    when :
    params.stats_alpha_enable
    
    script :
    """
    Rscript --vanilla ${baseDir}/lib/alpha_diversity.R phyloseq.rds ${params.stats.perc_abund_threshold} ${params.stats.distance} alpha_div_plots_${params.stats.alpha_div_group}.svg barplot_relabund_phylum_${params.stats.alpha_div_group}.svg barplot_relabund_family_${params.stats.alpha_div_group}.svg barplot_relabund_genus_${params.stats.alpha_div_group}.svg heatmap_class.svg heatmap_family.svg heatmap_genus.svg ${params.stats.alpha_div_group} > stats_alpha_diversity.log 2>&1
    cp ${baseDir}/lib/alpha_diversity.R completecmd >> stats_alpha_diversity.log 2>&1
    """
}

process stats_beta {

    beforeScript "${params.r_stats_env}"

    publishDir "${params.outdir}/${params.report_dirname}/R/SCRIPT", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_beta -> "${task.process}.R" }
    publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_non_normalized", mode: 'copy', pattern : '*.svg'
    publishDir "${params.outdir}/${params.report_dirname}/R/DATA", mode: 'copy', pattern : '*.tsv'

    input :
        file phyloseq_rds from phyloseq_rds_beta
        file metadata from metadata_beta
 
    output :
        file 'completecmd' into complete_cmd_beta
        file 'samples_ordination_plot_*.svg' into samples_ordination_plot

    //Run only if process is activated in params.config file
    when :
    params.stats_beta_enable

 
    script:
    """
    Rscript --vanilla ${baseDir}/lib/beta_diversity.R ${phyloseq_rds} ${params.stats.distance} ${params.stats.beta_div_criteria} samples_ordination_plot_${params.stats.beta_div_criteria}.svg ${metadata} $workflow.projectDir > stats_beta_diversity.log 2>&1
    cp ${baseDir}/lib/beta_diversity.R completecmd >> stats_beta_diversity.log 2>&1
    """
}

process stats_beta_rarefied {

    beforeScript "${params.r_stats_env}"
    publishDir "${params.outdir}/${params.report_dirname}/R/SCRIPT", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_beta_rarefied -> "${task.process}.R" }
    publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_rarefied", mode: 'copy', pattern : '*.svg'
    publishDir "${params.outdir}/${params.report_dirname}/R/DATA", mode: 'copy', pattern : '*.tsv'

    input :
        file phyloseq_rds from phyloseq_rds_beta_rarefied
        file metadata from metadata_beta_rarefied
     
    output :
        file 'completecmd' into complete_cmd_beta_rarefied
        file 'Final_rarefied_ASV_table_with_taxonomy.tsv' into final_rarefied_ASV_table_with_taxonomy
        file 'samples_ordination_plot_rarefied_*.svg' into samples_ordination_plot_rarefied

    //Run only if process is activated in params.config file
    when :
    params.stats_beta_enable

    script:
    """
    Rscript --vanilla ${baseDir}/lib/beta_diversity_rarefied.R ${phyloseq_rds} Final_rarefied_ASV_table_with_taxonomy.tsv ${params.stats.distance} ${params.stats.beta_div_criteria} samples_ordination_plot_rarefied_${params.stats.beta_div_criteria}.svg ${metadata} $workflow.projectDir > stats_beta_diversity_rarefied.log 2>&1
    cp ${baseDir}/lib/beta_diversity_rarefied.R completecmd >> stats_beta_diversity_rarefied.log 2>&1
    """
}

process stats_beta_deseq2 {

    beforeScript "${params.r_stats_env}"
    publishDir "${params.outdir}/${params.report_dirname}/R/SCRIPT", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_beta_deseq2 -> "${task.process}.R" }
    publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_deseq2", mode: 'copy', pattern : '*.svg'
    publishDir "${params.outdir}/${params.report_dirname}/R/DATA", mode: 'copy', pattern : '*.tsv'

    input :
        file phyloseq_rds from phyloseq_rds_beta_deseq2
        file metadata from metadata_beta_deseq2

    output :
        file 'completecmd' into complete_cmd_beta_deseq2
        file 'Final_deseq2_ASV_table_with_taxonomy.tsv' into final_deseq2_ASV_table_with_taxonomy
        file 'samples_ordination_plot_deseq2_*.svg' into samples_ordination_plot_deseq2

    //Run only if process is activated in params.config file
    when :
    params.stats_beta_enable

    script:
    """
    Rscript --vanilla ${baseDir}/lib/beta_diversity_deseq2.R ${phyloseq_rds} Final_deseq2_ASV_table_with_taxonomy.tsv ${params.stats.distance} ${params.stats.beta_div_criteria} samples_ordination_plot_deseq2_${params.stats.beta_div_criteria}.svg ${metadata} $workflow.projectDir > stats_beta_diversity_deseq2.log 2>&1
    cp ${baseDir}/lib/beta_diversity_deseq2.R completecmd >> stats_beta_diversity_deseq2.log 2>&1
    """
}

process stats_beta_css {

    beforeScript "${params.r_stats_env}"
    publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_beta_css -> "cmd/${task.process}.R" }
    publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_css", mode: 'copy', pattern : '*.svg'
    publishDir "${params.outdir}/${params.report_dirname}/R/DATA", mode: 'copy', pattern : '*.tsv'

    input :
        file phyloseq_rds from phyloseq_rds_beta_css
        file metadata from metadata_beta_css

    output :
        file 'completecmd' into complete_cmd_beta_css
        file 'Final_css_ASV_table_with_taxonomy.tsv' into final_css_ASV_table_with_taxonomy
        file 'samples_ordination_plot_css_*.svg' into samples_ordination_plot_css

    //Run only if process is activated in params.config file
    when :
    params.stats_beta_enable

    script:
    """
    Rscript --vanilla ${baseDir}/lib/beta_diversity_css.R ${phyloseq_rds} Final_css_ASV_table_with_taxonomy.tsv ${params.stats.distance} ${params.stats.beta_div_criteria} samples_ordination_plot_css_${params.stats.beta_div_criteria}.svg ${metadata} $workflow.projectDir > stats_beta_diversity_css.log 2>&1
    cp ${baseDir}/lib/beta_diversity_css.R completecmd >> stats_beta_diversity_css.log 2>&1
    """
} 
