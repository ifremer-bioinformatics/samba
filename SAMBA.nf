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
println "Data type : $params.data_type"

Channel.fromPath(params.inmanifest, checkIfExists:true).into { manifest ; manifest4integrity }
Channel.fromPath(params.inmetadata, checkIfExists:true).into { metadata; metadata_dbotu3 ; metadata4stats ; metadata4integrity ; metadata4picrust2 }

//Copy params.config file to output directory for each run
paramsfile = file('config/params.config')
paramsfile.copyTo("$params.outdir/config/params.config")

// IF NOT STATS ONLY, PERFORM QIIME STEPS
if(!params.stats_only) {
    /* Check data integrity */
    
    process data_integrity {
        publishDir "${params.outdir}/${params.data_integrity_dirname}", mode: 'copy', pattern: 'data_integrity.csv'
        publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: 'data_integrity.csv'

    input :
        file manifest from manifest4integrity
        file metadata from metadata4integrity

    output :
        file 'verifications.ok' optional true into ckeck_ok
        file 'verifications.bad' optional true into check_bad
        file 'data_integrity.csv' optional true into data_integrity_csv
    
    //Run only if process is activated in params.config file
    when :
        params.data_integrity_enable

    script :
    """
    ${baseDir}/lib/data_integrity.sh ${manifest} ${metadata} ${params.data_integrity.primerF} ${params.data_integrity.primerR} data_integrity.csv verifications.ok verifications.bad ${params.data_integrity.barcode_column_name} ${params.data_integrity.sampleid_column_name} ${params.data_integrity.R1_single_files_column_name} ${params.data_integrity.R1_files_column_name} ${params.data_integrity.R2_files_column_name} ${params.data_integrity.barcode_filter} ${params.data_integrity.primer_filter} ${params.data_type} &> data_integrity.log 2>&1
    if test -f "verifications.bad"; then
        if test -f "data_integrity.csv"; then
            echo "Data integrity process not satisfied, check ${params.outdir}/${params.data_integrity_dirname}/data_integrity.csv file"
            mkdir -p ${params.outdir}/${params.data_integrity_dirname}
            cp data_integrity.csv ${params.outdir}/${params.data_integrity_dirname}/.
        else
            echo "Data integrity process not satisfied, check data_integrity.log file in process working directory"
        fi
        exit 1
     fi

    """
   }

    /* Import metabarcode data */
    
    process q2_import {

        beforeScript "${params.load_conda}"
        conda "${params.qiime_env}" 

        publishDir "${params.outdir}/${params.import_dirname}", mode: 'copy', pattern: 'data.qz*'
        publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: '*_output'
        publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_import -> "cmd/${task.process}_complete.sh" }
    
        input : 
            file q2_manifest from manifest
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
        ${baseDir}/lib/q2_import.sh ${params.data_type} ${q2_manifest} data.qza data.qzv import_output completecmd &> q2_import.log 2>&1
        """
    }
    
    /* Trim metabarcode data with cutadapt */
    process q2_cutadapt {
    
        beforeScript "${params.load_conda}"
        conda "${params.qiime_env}" 

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
        ${baseDir}/lib/q2_cutadapt.sh ${params.data_type} ${task.cpus} ${imported_data} ${params.cutadapt.primerF} ${params.cutadapt.primerR} ${params.cutadapt.errorRate} ${params.cutadapt.overlap} data_trimmed.qza data_trimmed.qzv trimmed_output completecmd &> q2_cutadapt.log 2>&1
        """
    }

    /* Run dada2 */
    process q2_dada2 {
    
        beforeScript "${params.load_conda}"
        conda "${params.qiime_env}" 

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
        ${baseDir}/lib/q2_dada2.sh ${params.data_type} ${trimmed_data} ${metadata} rep_seqs.qza rep_seqs.qzv table.qza table.qzv stats.qza stats.qzv dada2_output ${params.dada2.trim5F} ${params.dada2.trim5R} ${params.dada2.trunclenF} ${params.dada2.trunclenR} ${params.dada2.maxee_f} ${params.dada2.maxee_r} ${params.dada2.minqual} ${params.dada2.chimeras} ${task.cpus} completecmd &> q2_dada2.log 2>&1
        """
    }
    
    data_table.set { table_picrust2 }
    data_repseqs.into { dada2_seqs_dbotu3 ; seqs_taxo ; seqs_phylo }
    dada2_summary.set { summary }

    if(params.dbotu3_enable){

    /* Run dbotu3 */
    process q2_dbotu3 {

        beforeScript "${params.load_conda}"
        conda "${params.dbotu_env}"

        publishDir "${params.outdir}/${params.dbotu3_dirname}", mode: 'copy', pattern: '*.qz*'
        publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: '*_output'
        publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_dbotu3 -> "cmd/${task.process}_complete.sh" }

        input :
            file table from data_table
            file seqs from dada2_seqs_dbotu3
            file metadata_dbotu3 from metadata_dbotu3

        output :
            file 'dbotu3_details.txt' into dbotu3_details
            file 'dbotu3_seqs.qza' into dbotu3_seqs
            file 'dbotu3_seqs.qzv' into dbotu3_seqs_visu
            file 'dbotu3_table.qza' into dbotu3_table
            file 'dbotu3_table.qzv' into dbotu3_table_visu
            file 'dbotu3_output' into dbotu3_summary
            file 'completecmd' into complete_cmd_dbotu3

        script :
        """
        ${baseDir}/lib/q2_dbotu3.sh ${table} ${seqs} ${metadata_dbotu3} dbotu3_details.txt dbotu3_seqs.qza dbotu3_seqs.qzv dbotu3_table.qza dbotu3_table.qzv dbotu3_output ${params.dbotu3.genet_crit} ${params.dbotu3.abund_crit} ${params.dbotu3.pval_crit} completecmd &> q2_dbotu3.log 2>&1
        """
    }

    dbotu3_table.set { table_picrust2 }
    dbotu3_seqs.into { seqs_taxo ; seqs_phylo ; seqs_picrust2 }
    dbotu3_summary.set { summary }
    }

    /* Run taxonomy assignment */

    process q2_taxonomy {
    
        beforeScript "${params.load_conda}"
        conda "${params.qiime_env}" 

        publishDir "${params.outdir}/${params.taxo_dirname}", mode: 'copy', pattern: '*.qz*'
        publishDir "${params.outdir}/${params.taxo_dirname}", mode: 'copy', pattern: '*.tsv*'
        publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: '*_output'
        publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: 'ASV_table*'
        publishDir "${params.outdir}/${params.report_dirname}/taxo_output/", mode: 'copy', pattern: 'ASV_taxonomy.tsv'
        publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_taxo -> "cmd/${task.process}_complete.sh" }
    
        input :
            file repseqs_taxo from seqs_taxo
            file summary from summary

        output :
            file 'taxonomy.qza' into data_taxonomy
            file 'taxonomy.qzv' into visu_taxonomy
            file 'ASV_taxonomy.tsv' into taxonomy_tsv
            file 'taxo_output' into taxo_summary
            file 'ASV_table_with_taxonomy.biom' into biom
            file 'ASV_table_with_taxonomy.tsv' into biom_tsv
            file 'taxonomic_database.qza' optional true into trained_database
            file 'db_seqs_amplicons.qza' optional true into db_seqs_filtered
            file 'completecmd' into complete_cmd_taxo
    
        //Run only if process is activated in params.config file
        when :
        params.taxo_enable
    
        script :
        """
        ${baseDir}/lib/q2_taxo.sh ${task.cpus} ${params.taxo.db_seqs} ${params.taxo.db_tax} ${params.taxo.database} ${params.taxo.extract_db} ${params.cutadapt.primerF} ${params.cutadapt.primerR} ${params.taxo.confidence} ${repseqs_taxo} taxonomy.qza taxonomy.qzv taxo_output ASV_taxonomy.tsv ${summary} ASV_table_with_taxonomy.biom ASV_table_with_taxonomy.tsv taxonomic_database.qza db_seqs_amplicons.qza completecmd &> q2_taxo.log 2>&1
        """ 
    }

    biom_tsv.set { tsv }

    if(params.microDecon_enable) {

        /* Run sample decontamination using MicroDecon */
    
        process microDecon_step1 {
    
            beforeScript "${params.load_conda}"
            conda "${params.microdecon_env}"
    
            publishDir "${params.outdir}/${params.microDecon_dirname}", mode: 'copy', pattern: 'decontaminated_ASV_table.tsv'
            publishDir "${params.outdir}/${params.microDecon_dirname}", mode: 'copy', pattern: 'abundance_removed.txt'
            publishDir "${params.outdir}/${params.microDecon_dirname}", mode: 'copy', pattern: 'ASV_removed.txt'
            publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: 'decontaminated_ASV_table.tsv'
            publishDir "${params.outdir}/${params.report_dirname}/microDecon", mode: 'copy', pattern: 'decontaminated_ASV_table.tsv'
            publishDir "${params.outdir}/${params.report_dirname}/microDecon", mode: 'copy', pattern: 'abundance_removed.txt'
            publishDir "${params.outdir}/${params.report_dirname}/microDecon", mode: 'copy', pattern: 'ASV_removed.txt'
            publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_microDecon -> "cmd/${task.process}_complete.sh" }
    
            input :
                file microDecon_table from biom_tsv
    
            output :
                file 'decontaminated_ASV_table.tsv' into decontam_table
                file 'abundance_removed.txt' into abund_removed
                file 'ASV_removed.txt' into ASV_removed
                file 'completecmd' into complete_cmd_microDecon
    
            shell :
            """ 
            sed '1d' ${microDecon_table} > microDecon_table
            sed -i 's/#OTU ID/ASV_ID/g' microDecon_table
            ${baseDir}/lib/microDecon.R microDecon_table ${params.microDecon.control_list} ${params.microDecon.nb_controls} ${params.microDecon.nb_samples} decontaminated_ASV_table.tsv abundance_removed.txt ASV_removed.txt &> microDecon.log 2>&1
            cp ${baseDir}/lib/microDecon.R completecmd &>> microDecon.log 2>&1
            
            """
        }

        decontam_table.into { decontam_table_step2 ; decontam_table_step3 ; tsv }

        process microDecon_step2 {

            beforeScript "${params.load_conda}"
            conda "${params.qiime_env}"

            publishDir "${params.outdir}/${params.microDecon_dirname}", mode: 'copy', pattern: 'decontaminated_ASV_table.qza'
            publishDir "${params.outdir}/${params.report_dirname}/microDecon", mode: 'copy', pattern: 'decontaminated_ASV_table.qza'

            input :
                file table4microDecon from decontam_table_step2

            output :
                file 'decontaminated_ASV_table.qza' into decontam_table_qza

            shell :
            """
            biom convert -i ${table4microDecon} -o decontaminated_ASV_table.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy
            qiime tools import --input-path decontaminated_ASV_table.biom --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path decontaminated_ASV_table.qza
            """
        }

        decontam_table_qza.set { table_picrust2 }
    
        /* Extract non-contaminated ASV ID */
    
        process microDecon_step3 {
    
            beforeScript "${params.load_conda}"
            conda "${params.seqtk_env}"
    
            publishDir "${params.outdir}/${params.microDecon_dirname}", mode: 'copy', pattern: 'decontaminated_ASV_ID.txt'
            publishDir "${params.outdir}/${params.microDecon_dirname}", mode: 'copy', pattern: 'decontaminated_ASV.fasta'
            publishDir "${params.outdir}/${params.report_dirname}/microDecon", mode: 'copy', pattern: 'decontaminated_ASV.fasta'
    
            input :
                file decontam_table from decontam_table_step3
                file dada2_summary from dada2_summary
    
            output :
                file 'decontaminated_ASV_ID.txt' into decontam_ASV_ID
                file 'decontaminated_ASV.fasta' into decontam_ASV_fasta
    
            shell :
            """
            cut -d \$'\t' -f1 ${decontam_table} | sed '1d' > decontaminated_ASV_ID.txt
            seqtk subseq ${dada2_summary}/sequences.fasta decontaminated_ASV_ID.txt > decontaminated_ASV.fasta
            """
        }

        decontam_ASV_fasta.set { seqs_phylo }
    
        /* Run phylogeny from decontaminated ASV sequences */
    
        process microDecon_step4 {
    
            beforeScript "${params.load_conda}"
            conda "${params.qiime_env}"
            
            publishDir "${params.outdir}/${params.phylogeny_dirname}", mode: 'copy', pattern: '*.qza'
            publishDir "${params.outdir}/${params.phylogeny_dirname}", mode: 'copy', pattern: '*.txt'
            publishDir "${params.outdir}/${params.phylogeny_dirname}", mode: 'copy', pattern: 'tree_export_dir'
            publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: 'tree_export_dir'
            publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_phylo -> "cmd/${task.process}_complete.sh" }
    
    
            input :
                file decontam_ASV_fasta from seqs_phylo
    
            output :
                file 'decontam_seqs.qza' into decontam_seqs_qza
                file 'aligned_repseq.qza' into aligned_repseq
                file 'masked-aligned_repseq.qza' into masked_aligned
                file 'tree.qza' into tree
                file 'tree.log' into tree_bestmodel_log
                file 'rooted_tree.qza' into rooted_tree
                file 'tree_export_dir' into tree_export_dir
                file 'tree_export.log' into tree_export_log
                file 'tree.nwk' into newick
                file 'completecmd' into complete_cmd_phylogeny
    
            shell :
            """
            qiime tools import --input-path ${decontam_ASV_fasta} --output-path decontam_seqs.qza --type 'FeatureData[Sequence]'
            ${baseDir}/lib/q2_phylogeny.sh decontam_seqs.qza aligned_repseq.qza masked-aligned_repseq.qza tree.qza tree.log rooted_tree.qza tree_export_dir tree_export.log completecmd &> q2_phylogeny.log 2>&1
            cp tree_export_dir/tree.nwk tree.nwk &>> q2_phylogeny.log 2>&1
            
            """
        }
    decontam_seqs_qza.set { seqs_picrust2 } 
    }    
    else {
    
        /* Run phylogeny construction */
    
        process q2_phylogeny {
    
            beforeScript "${params.load_conda}"
            conda "${params.qiime_env}"
    
            publishDir "${params.outdir}/${params.phylogeny_dirname}", mode: 'copy', pattern: '*.qza'
            publishDir "${params.outdir}/${params.phylogeny_dirname}", mode: 'copy', pattern: '*.txt'
            publishDir "${params.outdir}/${params.phylogeny_dirname}", mode: 'copy', pattern: 'tree_export_dir'
            publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: 'tree_export_dir'
            publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_phylo -> "cmd/${task.process}_complete.sh" }
    
            input :
                file repseqs_phylo from seqs_phylo
    
            output :
                file 'aligned_repseq.qza' into aligned_repseq
                file 'masked-aligned_repseq.qza' into masked_aligned
                file 'tree.qza' into tree
                file 'tree.log' into tree_bestmodel_log
                file 'rooted_tree.qza' into rooted_tree
                file 'tree_export_dir' into tree_export_dir
                file 'tree_export.log' into tree_export_log
                file 'tree.nwk' into newick
                file 'completecmd' into complete_cmd_phylogeny
    
            //Run only if process is activated in params.config file
            when :
            params.phylogeny_enable
    
            script :
            """
            ${baseDir}/lib/q2_phylogeny.sh ${repseqs_phylo} aligned_repseq.qza masked-aligned_repseq.qza tree.qza tree.log rooted_tree.qza tree_export_dir tree_export.log completecmd &> q2_phylogeny.log 2>&1
            cp tree_export_dir/tree.nwk tree.nwk &>> q2_phylogeny.log 2>&1
            """
        }
    }

    /* Run functional predictions */

    process q2_picrust2_analysis {

        beforeScript "${params.load_conda}"
        conda "${params.picrust_env}"

        publishDir "${params.outdir}/${params.picrust2_dirname}", mode: 'copy', pattern: 'q2-picrust2_output/*'
        publishDir "${params.outdir}/${params.picrust2_dirname}", mode: 'copy', pattern: 'q2-picrust2_output/*_exported/*.tsv'
        publishDir "${params.outdir}/${params.report_dirname}/picrust2_output", mode: 'copy', pattern: 'q2-picrust2_output/*'
        publishDir "${params.outdir}/${params.report_dirname}/picrust2_output", mode: 'copy', pattern: 'q2-picrust2_output/*_exported/*.tsv'
        publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'complete_picrust2_cmd', saveAs : { complete_picrust2_cmd -> "cmd/${task.process}_complete.sh" }

        input :
            file seqs_picrust2 from seqs_picrust2
            file table_picrust2 from table_picrust2

        output :
            file 'q2-picrust2_output/ec_metagenome.qza' into EC_predictions
            file 'q2-picrust2_output/ec_metagenome.qzv' into EC_predictions_visu
            file 'q2-picrust2_output/ec_metagenome_exported/ec_metagenome_predictions*.tsv' into EC_predictions_tsv
            file 'q2-picrust2_output/ko_metagenome.qza' into KO_predictions
            file 'q2-picrust2_output/ko_metagenome.qzv' into KO_predictions_visu
            file 'q2-picrust2_output/ko_metagenome_exported/ko_metagenome_predictions*.tsv' into KO_predictions_tsv
            file 'q2-picrust2_output/pathway_abundance.qza' into pathway_predictions
            file 'q2-picrust2_output/pathway_abundance_visu' into pathway_predictions_visu
            file 'q2-picrust2_output/pathway_abundance_exported/pathway_abundance_predictions*.tsv' into pathway_predictions_tsv
            file 'complete_picrust2_cmd' into complete_picrust2_cmd

        //Run only if process is activated in params.config file
        when :
        params.picrust2_enable

        script :
        """
        ${baseDir}/lib/q2_picrust2.sh ${table_picrust2} ${seqs_picrust2} q2-picrust2_output ${task.cpus} ${params.picrust2.method} ${params.picrust2.nsti} complete_picrust2_cmd &> q2_picrust2.log 2>&1
        """
    }
    
    /* Statistical analysis of functional predictions  */

    process q2_picrust2_stats {

    beforeScript "${params.load_conda}"
    conda "${params.r_stats_env}"

    publishDir "${params.outdir}/${params.report_dirname}/picrust2_output", mode: 'copy', pattern: '*functional_predictions_NMDS*'
    publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'complete_picrust2_stats_cmd', saveAs : { complete_picrust2_stats_cmd -> "cmd/${task.process}_complete.sh" }

    input :
        file ec_metagenome from EC_predictions_tsv
        file ko_metagenome from KO_predictions_tsv
        file metacyc_predictions_ from pathway_predictions_tsv
        file metadata4picrust2 from metadata4picrust2
    output :
        file '*functional_predictions_NMDS*' into functional_pred_NMDS
        file 'complete_picrust2_stats_cmd' into complete_picrust2_stats_cmd

    //Run only if process is activated in params.config file
    when :
    params.picrust2_enable

    script :
    """
    ${baseDir}/lib/functional_predictions.R ec_metagenome_predictions_with-descriptions.tsv ko_metagenome_predictions_with-descriptions.tsv pathway_abundance_predictions_with-descriptions.tsv ${metadata4picrust2} ${params.stats.beta_div_criteria} functional_predictions_NMDS ${params.microDecon_enable} ${params.microDecon.control_list} &> picrust2_stats.log 2>&1
    cp ${baseDir}/lib/functional_predictions.R complete_picrust2_stats_cmd &>> picrust2_stats.log 2>&1
    """
    }
}

if(params.stats_only){
    
    //IF OTU TABLE ALREADY CREATED AND STAT ONLY STEPS NEEDED
    
    //Set biom_tsv path in params.conf
    Channel.fromPath(params.newick, checkIfExists:true).set { newick }
    Channel.fromPath(params.inasv_table, checkIfExists:true).set { tsv }
    println "Input ASV table used for statistics steps : $params.inasv_table"
}


process prepare_data_for_stats {

    beforeScript "${params.load_conda}"
    conda "${params.r_stats_env}"

    publishDir "${params.outdir}/${params.report_dirname}/R/DATA", mode: 'copy', pattern : '*.tsv'
    publishDir "${params.outdir}/${params.report_dirname}/R/DATA", mode: 'copy', pattern : '*.rds'
    
    input :
        file metadata from metadata4stats
        file biom_tsv from tsv
        file newick from newick

    output :
        file 'ASV_table_with_taxo_for_stats.tsv' into biom_tsv_stats
        file 'metadata_stats.tsv' into metadata_stats
        file 'phyloseq.rds' into phyloseq_rds
 
    when :
    params.prepare_data_for_stats_enable

    script :
    """
    ${baseDir}/lib/prepare_data_for_stats.sh ${metadata} ${biom_tsv} ASV_table_with_taxo_for_stats.tsv metadata_stats.tsv ${params.microDecon_enable} &> stats_prepare_data.log 2&>1
    Rscript --vanilla ${baseDir}/lib/create_phyloseq_obj.R phyloseq.rds ASV_table_with_taxo_for_stats.tsv metadata_stats.tsv ${params.microDecon_enable} ${params.microDecon.control_list} ${newick} &>> stats_prepare_data.log 2&>1 
    """
}

//Duplicate channels needed in several processes
metadata_stats.into { metadata_beta ; metadata_beta_rarefied ; metadata_beta_deseq2 ; metadata_beta_css }
phyloseq_rds.into { phyloseq_rds_alpha ; phyloseq_rds_beta ; phyloseq_rds_beta_rarefied ; phyloseq_rds_beta_deseq2 ; phyloseq_rds_beta_css }

process stats_alpha {

    beforeScript "${params.load_conda}"
    conda "${params.r_stats_env}"

    publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/alpha_diversity", mode: 'copy', pattern : 'index_significance_tests.txt'
    publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/alpha_diversity/diversity_index", mode: 'copy', pattern : 'alpha_div_plots*'
    publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/alpha_diversity/diversity_barplots", mode: 'copy', pattern : 'barplot_*'
    publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/alpha_diversity", mode: 'copy', pattern : 'rarefaction_curve*'
    
    input :
        file phyloseq_rds from phyloseq_rds_alpha

    output :
        file 'alpha_div_plots*' into alpha_div_plots
        file 'index_significance_tests.txt' into index_significance_tests
        file 'barplot_phylum*' into barplot_phylum
        file 'barplot_class*' into barplot_class
        file 'barplot_order*' into barplot_order
        file 'barplot_family*' into barplot_family
        file 'barplot_genus*' into barplot_genus
        file 'rarefaction_curve*' into rarefaction_curve

    //Run only if process is activated in params.config file
    when :
    params.stats_alpha_enable
    
    script :
    """
    Rscript --vanilla ${baseDir}/lib/alpha_diversity.R phyloseq.rds ${params.stats.distance} alpha_div_plots ${params.stats.kingdom} ${params.stats.nbtax} barplot_phylum barplot_class barplot_order barplot_family barplot_genus ${params.stats.alpha_div_group} index_significance_tests.txt $workflow.projectDir rarefaction_curve &> stats_alpha_diversity.log 2>&1
    """
}

process stats_beta {

    beforeScript "${params.load_conda}"
    conda "${params.r_stats_env}"

    publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_non_normalized/NMDS", mode: 'copy', pattern : 'NMDS*'
    publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_non_normalized/PCoA", mode: 'copy', pattern : 'PCoA*'
    publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_non_normalized/Hierarchical_Clustering", mode: 'copy', pattern : 'hclustering*'
    publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_non_normalized", mode:'copy', pattern : 'variance_significance_tests_*'
    publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_non_normalized", mode:'copy', pattern : 'pie_ExpVar_*'
    publishDir "${params.outdir}/${params.report_dirname}/R/DATA", mode: 'copy', pattern : '*.tsv'

    input :
        file phyloseq_rds from phyloseq_rds_beta
        file metadata from metadata_beta
 
    output :
        file 'NMDS_*' into NMDS
        file 'PCoA_*' into PCoA
        file 'hclustering_*' into hclustering
        file 'variance_significance_tests_*' into variance_significance_tests
        file 'pie_ExpVar_*' into pie_ExpVar

    //Run only if process is activated in params.config file
    when :
    params.stats_beta_enable

 
    script:
    """
    Rscript --vanilla ${baseDir}/lib/beta_diversity.R ${phyloseq_rds} ${params.stats.beta_div_criteria} ${metadata} $workflow.projectDir NMDS_ PCoA_ ${params.stats.hc_method} hclustering_ variance_significance_tests_ pie_ExpVar_ &> stats_beta_diversity.log 2>&1
    """
}

process stats_beta_rarefied {

    beforeScript "${params.load_conda}"
    conda "${params.r_stats_env}"

    publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_rarefied/NMDS", mode: 'copy', pattern : 'NMDS*'
    publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_rarefied/PCoA", mode: 'copy', pattern : 'PCoA*'
    publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_rarefied/Hierarchical_Clustering", mode: 'copy', pattern : 'hclustering*'
    publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_rarefied", mode: 'copy', pattern : 'variance_significance_tests_rarefied_*'
    publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_rarefied", mode:'copy', pattern : 'pie_ExpVar_rarefied_*'
    publishDir "${params.outdir}/${params.report_dirname}/R/DATA", mode: 'copy', pattern : '*.tsv'

    input :
        file phyloseq_rds from phyloseq_rds_beta_rarefied
        file metadata from metadata_beta_rarefied
     
    output :
        file 'Final_rarefied_ASV_table_with_taxonomy.tsv' into final_rarefied_ASV_table_with_taxonomy
        file 'NMDS_rarefied_*' into NMDS_rarefied
        file 'PCoA_rarefied_*' into PCoA_rarefied
        file 'hclustering_rarefied_*' into hclustering_rarefied
        file 'variance_significance_tests_rarefied_*' into variance_significance_tests_rarefied
        file 'pie_ExpVar_rarefied_*' into pie_ExpVar_rarefied

    //Run only if process is activated in params.config file
    when :
    params.stats_beta_enable

    script:
    """
    Rscript --vanilla ${baseDir}/lib/beta_diversity_rarefied.R ${phyloseq_rds} Final_rarefied_ASV_table_with_taxonomy.tsv ${params.stats.beta_div_criteria} ${metadata} $workflow.projectDir NMDS_rarefied_ PCoA_rarefied_ ${params.stats.hc_method} hclustering_rarefied_ variance_significance_tests_rarefied_ pie_ExpVar_rarefied_ &> stats_beta_diversity_rarefied.log 2>&1
    """
}

process stats_beta_deseq2 {

    beforeScript "${params.load_conda}"
    conda "${params.r_stats_env}"

    publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_DESeq2/NMDS", mode: 'copy', pattern : 'NMDS*'
    publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_DESeq2/PCoA", mode: 'copy', pattern : 'PCoA*'
    publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_DESeq2/Hierarchical_Clustering", mode: 'copy', pattern : 'hclustering*'
    publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_DESeq2", mode: 'copy',pattern : 'variance_significance_tests_DESeq2_*'
    publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_DESeq2", mode: 'copy',pattern : 'pie_ExpVar_DESeq2_*'
    publishDir "${params.outdir}/${params.report_dirname}/R/DATA", mode: 'copy', pattern : '*.tsv'

    input :
        file phyloseq_rds from phyloseq_rds_beta_deseq2
        file metadata from metadata_beta_deseq2

    output :
        file 'Final_DESeq2_ASV_table_with_taxonomy.tsv' into final_deseq2_ASV_table_with_taxonomy
        file 'NMDS_DESeq2_*' into NMDS_deseq2
        file 'PCoA_DESeq2_*' into PCoA_deseq2
        file 'hclustering_DESeq2_*' into hclustering_deseq2
        file 'variance_significance_tests_DESeq2_*' into variance_significance_tests_DESeq2
        file 'pie_ExpVar_DESeq2_*' into pie_ExpVar_DESeq2
        
    //Run only if process is activated in params.config file
    when :
    params.stats_beta_enable

    script:
    """
    Rscript --vanilla ${baseDir}/lib/beta_diversity_deseq2.R ${phyloseq_rds} Final_DESeq2_ASV_table_with_taxonomy.tsv ${params.stats.beta_div_criteria} ${metadata} $workflow.projectDir NMDS_DESeq2_ PCoA_DESeq2_ ${params.stats.hc_method} hclustering_DESeq2_ variance_significance_tests_DESeq2_ pie_ExpVar_DESeq2_ &> stats_beta_diversity_deseq2.log 2>&1
    """
}

process stats_beta_css {

    beforeScript "${params.load_conda}"
    conda "${params.r_stats_env}"

    publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_CSS/NMDS", mode: 'copy', pattern : 'NMDS*'
    publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_CSS/PCoA", mode: 'copy', pattern : 'PCoA*'
    publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_CSS/Hierarchical_Clustering", mode: 'copy', pattern : 'hclustering*'
    publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_CSS", mode: 'copy', pattern : 'variance_significance_tests_CSS_*'
    publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_CSS", mode: 'copy', pattern : 'pie_ExpVar_CSS_*'
    publishDir "${params.outdir}/${params.report_dirname}/R/DATA", mode: 'copy', pattern : '*.tsv'

    input :
        file phyloseq_rds from phyloseq_rds_beta_css
        file metadata from metadata_beta_css

    output :
        file 'Final_CSS_ASV_table_with_taxonomy.tsv' into final_css_ASV_table_with_taxonomy
        file 'NMDS_CSS_*' into NMDS_css
        file 'PCoA_CSS_*' into PCoA_css
        file 'hclustering_CSS_*' into hclustering_css
        file 'variance_significance_tests_CSS_*' into variance_significance_tests_CSS
        file 'pie_ExpVar_CSS_*' into pie_ExpVar_CSS
        file 'end_analysis.ok' into end_ok


    //Run only if process is activated in params.config file
    when :
    params.stats_beta_enable

    script:
    """
    Rscript --vanilla ${baseDir}/lib/beta_diversity_css.R ${phyloseq_rds} Final_CSS_ASV_table_with_taxonomy.tsv ${params.stats.beta_div_criteria} ${metadata} $workflow.projectDir NMDS_CSS_ PCoA_CSS_ ${params.stats.hc_method} hclustering_CSS_ variance_significance_tests_CSS_ pie_ExpVar_CSS_ &> stats_beta_diversity_css.log 2>&1
    touch end_analysis.ok
    """
}

Channel.fromPath(params.reportHTML, checkIfExists:true).set { reportHTML }
Channel.fromPath(params.reportMD, checkIfExists:true).set { reportMD }

process report {

    publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'Report_*'

    input :
        file reportHTML from reportHTML
        file reportMD from reportMD
        file 'end_analysis.ok' from end_ok

    output :
        file 'Report_*' into Reports

    //Run only if process is activated in params.config file
    when :
    params.report_enable

    script:
    """
    cp ${reportHTML} Report_${params.projectName}.html
    cp ${reportMD} Report_${params.projectName}.md
    """
}
