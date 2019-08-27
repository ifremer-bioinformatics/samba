#!/usr/bin/env nextflow
println "Workflow for project : $params.projectName"
println "Workflow description : $workflow.manifest.description"
println "Workflow gitLab URL : $workflow.manifest.homePage"
println "Workflow authors : $workflow.manifest.author"
println "Project : $workflow.projectDir"
println "Cmd line: $workflow.commandLine"
println "Workflow working/temp directory : $workflow.workDir"
println "Workflow output/publish directory : $params.outdir"
println "Workflow configuration file : $workflow.configFiles"

Channel.fromPath(params.inmanifest, checkIfExists:true).set { manifest }
Channel.fromPath(params.inmetadata, checkIfExists:true).set { metadata }
Channel.fromPath(params.inmetadata, checkIfExists:true).set { metadata4stats }

Channel.fromPath(params.stats.Rstatsscript, checkIfExists:true).set { Rstatsscript }

/* Import metabarcode data */
process q2_import {

    beforeScript "${params.qiime_env}"
    publishDir "${params.outdir}/01_import", mode: 'copy', pattern: 'data.qz*'
    publishDir "${params.outdir}/00_report", mode: 'copy', pattern: '*_output'
    publishDir "${params.outdir}/00_report", mode: 'copy', pattern : '.command.sh', saveAs : { cmd_import -> "cmd/${task.process}.sh" }

    input : 
        file manifest from manifest

    output : 
        file 'data.qza' into imported_data
        file 'data.qzv' into imported_visu
        file 'import_output' into imported_summary
        file '.command.sh' into cmd_import

    script :
    """
    ${baseDir}/lib/q2_import.sh ${manifest} data.qza data.qzv import_output > q2_import.log 2>&1
    """
}

/* Trim metabarcode data with cutadapt */
process q2_cutadapt {

    beforeScript "${params.qiime_env}"
    publishDir "${params.outdir}/02_trimmed", mode: 'copy', pattern: 'data*.qz*'
    publishDir "${params.outdir}/00_report", mode: 'copy', pattern: '*_output'
    publishDir "${params.outdir}/00_report", mode: 'copy', pattern : '.command.sh', saveAs : { cmd_cutadapt -> "cmd/${task.process}.sh" }
    
    input : file imported_data from imported_data

    output :
        file 'data_trimmed.qza' into trimmed_data
        file 'data_trimmed.qzv' into trimmed_visu
        file 'trimmed_output' into trimmed_summary
        file '.command.sh' into cmd_cutadapt 

    //Run only if process is activated in params.config file
    when :
    params.cutadapt.enable

    script :
    """
    #Run cutadapt
    qiime cutadapt trim-paired \
    --verbose \
    --p-cores ${task.cpus} \
    --i-demultiplexed-sequences ${imported_data} \
    --p-front-f ${params.cutadapt.primerF} \
    --p-front-r ${params.cutadapt.primerR} \
    --p-error-rate ${params.cutadapt.errorRate} \
    --p-discard-untrimmed --p-match-read-wildcards \
    --p-overlap ${params.cutadapt.overlap} \
    --o-trimmed-sequences data_trimmed.qza > q2_cutadapt.log 
    #Summarize counts per sample for all samples
    qiime demux summarize \
    --verbose \
    --i-data data_trimmed.qza \
    --o-visualization data_trimmed.qzv >> q2_cutadapt.log 2>&1
    #Export html report
    qiime tools export \
    --input-path data_trimmed.qzv \
    --output-path trimmed_output >> q2_cutadapt.log 2>&1
    """
}

/* Run dada2 */
process q2_dada2 {

    beforeScript "${params.qiime_env}"
    publishDir "${params.outdir}/03_dada2", mode: 'copy', pattern: '*.qz*'
    publishDir "${params.outdir}/00_report", mode: 'copy', pattern: '*_output'
    publishDir "${params.outdir}/00_report", mode: 'copy', pattern : '.command.sh', saveAs : { cmd_dada2 -> "cmd/${task.process}.sh" }

    input : 
        file trimmed_data from trimmed_data
        file metadata from metadata        

    output :
        file 'rep_seqs.qza' into data_repseqs
        file 'table.qza' into data_table
        file 'stats.qza' into data_stats
        file 'rep_seqs.qzv' into visu_repseqs
        file 'table.qzv' into visu_table
        file 'stats.qzv' into visu_stats
        file 'dada2_output' into dada2_summary
        file '.command.sh' into cmd_dada2
    
    //Run only if process is activated in params.config file
    when :
    params.dada2.enable

    script :
    """
    #Run dada2 : denoises paired-end sequences, dereplicates them and filters chimeras
    qiime dada2 denoise-paired \
    --verbose \
    --i-demultiplexed-seqs ${trimmed_data} \
    --p-trim-left-f ${params.dada2.trim3F} \
    --p-trim-left-r ${params.dada2.trim3R} \
    --p-trunc-len-f ${params.dada2.trunclenF} \
    --p-trunc-len-r ${params.dada2.trunclenR} \
    --p-max-ee ${params.dada2.maxee} \
    --p-trunc-q ${params.dada2.minqual} \
    --p-chimera-method ${params.dada2.chimeras} \
    --p-n-threads ${task.cpus} \
    --o-representative-sequences rep_seqs.qza \
    --o-table table.qza \
    --o-denoising-stats stats.qza > q2_dada2.log 2>&1
    #Generate a tabular view of taxonomy metadata
    qiime metadata tabulate \
    --verbose \
    --m-input-file stats.qza \
    --o-visualization stats.qzv >> q2_dada2.log 2>&1
    #Generate visual and tabular summaries of a feature table
    qiime feature-table summarize \
    --verbose \
    --i-table table.qza \
    --o-visualization table.qzv \
    --m-sample-metadata-file ${metadata} >> q2_dada2.log 2>&1
    #Generate tabular view of feature identifier to sequence mapping
    qiime feature-table tabulate-seqs \
    --verbose \
    --i-data rep_seqs.qza \
    --o-visualization rep_seqs.qzv >> q2_dada2.log 2>&1
    #Export data to html
    qiime tools export \
    --input-path rep_seqs.qzv \
    --output-path dada2_output >> q2_dada2.log 2>&1
    qiime tools export \
    --input-path table.qzv \
    --output-path dada2_output >> q2_dada2.log 2>&1
    qiime tools export \
    --input-path stats.qzv \
    --output-path dada2_output >> q2_dada2.log 2>&1
    qiime tools export \
    --input-path table.qza \
    --output-path dada2_output >> q2_dada2.log 2>&1
    """
}

/* Run taxonomy assignment */
process q2_taxonomy {

    beforeScript "${params.qiime_env}"
    publishDir "${params.outdir}/04_taxonomy", mode: 'copy', pattern: '*.qz*'
    publishDir "${params.outdir}/04_taxonomy", mode: 'copy', pattern: '*.tsv*'
    publishDir "${params.outdir}/00_report", mode: 'copy', pattern: '*_output'
    publishDir "${params.outdir}/00_report", mode: 'copy', pattern: 'Final_ASV_table*'
    publishDir "${params.outdir}/00_report", mode: 'copy', pattern : '.command.sh', saveAs : { cmd_taxo -> "cmd/${task.process}.sh" }

    input :
        file data_repseqs from data_repseqs
        file dada2_summary from dada2_summary

    output :
        file 'taxonomy.qza' into data_taxonomy
        file 'taxonomy.qzv' into visu_taxonomy
        file 'ASV_taxonomy.tsv' into taxonomy_tsv
        file 'taxo_output' into taxo_summary
        file 'Final_ASV_table_with_taxonomy.biom' into biom
        file 'Final_ASV_table_with_taxonomy.tsv' into biom_tsv
        file '.command.sh' into cmd_taxo

    //Run only if process is activated in params.config file
    when :
    params.taxo.enable

    script :
    """
    #Run RDP Classifier for taxonomy assignment
    qiime feature-classifier classify-sklearn \
    --p-n-jobs ${task.cpus} \
    --p-confidence ${params.taxo.confidence} \
    --i-classifier ${params.taxo.database} \
    --i-reads ${data_repseqs} \
    --o-classification taxonomy.qza > q2_taxo.log 2>&1
    #Generate a tabular view of taxonomy metadata
    qiime metadata tabulate \
    --m-input-file taxonomy.qza \
    --o-visualization taxonomy.qzv >> q2_taxo.log 2>&1
    #Export data
    qiime tools export \
    --input-path taxonomy.qzv \
    --output-path taxo_output >> q2_taxo.log 2>&1
    #Rename tabular taxonomy file and modify header
    mv taxo_output/metadata.tsv ASV_taxonomy.tsv >> q2_taxo.log 2>&1
    sed -i '1,2d' ASV_taxonomy.tsv >> q2_taxo.log 2>&1
    sed -i '1 i\\#OTUID\ttaxonomy\tconfidence' ASV_taxonomy.tsv >> q2_taxo.log 2>&1
    #Add taxonomy to count table (biom format)
    biom add-metadata \
    -i ${dada2_summary}/feature-table.biom \
    --observation-metadata-fp ASV_taxonomy.tsv \
    -o Final_ASV_table_with_taxonomy.biom \
    --sc-separated taxonomy >> q2_taxo.log 2>&1
    #Convert biom table to tabular
    biom convert \
    -i Final_ASV_table_with_taxonomy.biom \
    -o Final_ASV_table_with_taxonomy.tsv \
    --to-tsv \
    --header-key taxonomy >> q2_taxo.log 2>&1
    """ 
}

/*
process statisticalanalysis {

    beforeScript "${params.r_stats_env}"
    publishDir "${params.outdir}/05_statistical_analysis", mode: 'copy'
    publishDir "${params.outdir}/00_report", mode: 'copy', pattern : '.command.sh', saveAs : { cmd_taxo -> "cmd/${task.process}.sh" }
    
    input :
        file Rstatsscript from Rstatsscript
        file biom_tsv from biom_tsv
        file metadata_stats from metadata4stats

    output :
        //file 'barplot_relabund_phylum.svg' into barplot_relabund_phylym
        //file 'barplot_relabund_phylum_rarefied.svg' into barplot_relabund_phylum_rarefied

    //Run only if process is activated in params.config file
    when :
    params.stats.enable
    
    script :
    """
    mkdir -p ${params.stats.data} ${params.stats.script} ${params.stats.figures}
    cp ${metadata_stats} ${params.stats.data}
    sed -i 's/#SampleID/SampleID/g' ${params.stats.data}/${metadata_stats}
    cp ${biom_tsv} ${params.stats.data}
    sed -i '1d' ${params.stats.data}/${biom_tsv}
    sed -i 's/#OTU ID/ASV_ID/g' ${params.stats.data}/${biom_tsv}
    Rscript --vanilla ${Rstatsscript} "toto" ${params.projectName} ${params.stats.perc_abund_threshold} ${params.stats.distance} ${params.stats.column_sample_replicat} ${params.stats.exp_var_samples}
    """
}
*/
