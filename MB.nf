#!/usr/bin/env nextflow
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
Channel.fromPath(params.stats.Rstatsscript, checkIfExists:true).set { Rstatsscript }

/* Import metabarcode data */
process q2_import {

    beforeScript '. /appli/bioinfo/qiime/2019.04/env.sh'
    publishDir "${params.outdir}/01_import", mode: 'copy', pattern: 'data.qz*'
    publishDir "${params.outdir}/00_report", mode: 'copy', pattern: '*_output'

    input : file manifest from manifest

    output : 
        file 'data.qza' into imported_data
        file 'data.qzv' into imported_visu
        file 'import_output' into imported_summary

    shell :
    """
    qiime tools import \
    --input-path ${manifest} \
    --output-path data.qza \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-format PairedEndFastqManifestPhred33V2 > q2_import.log 2>&1
    qiime demux summarize \
    --verbose \
    --i-data data.qza \
    --o-visualization data.qzv >> q2_import.log 2>&1
    qiime tools export \
    --input-path data.qzv \
    --output-path import_output >> q2_import.log 2>&1
    """
}

/* Trim metabarcode data with cutadapt */
process q2_cutadapt {

    beforeScript '. /appli/bioinfo/qiime/2019.04/env.sh'
    publishDir "${params.outdir}/02_trimmed", mode: 'copy', pattern: 'data*.qz*'
    publishDir "${params.outdir}/00_report", mode: 'copy', pattern: '*_output'
    
    input : file imported_data from imported_data

    output :
        file 'data_trimmed.qza' into trimmed_data
        file 'data_trimmed.qzv' into trimmed_visu
        file 'trimmed_output' into trimmed_summary

    shell :
    """
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
    qiime demux summarize \
    --verbose \
    --i-data data_trimmed.qza \
    --o-visualization data_trimmed.qzv >> q2_cutadapt.log 2>&1
    qiime tools export \
    --input-path data_trimmed.qzv \
    --output-path trimmed_output >> q2_cutadapt.log 2>&1
    """
}

/* Run dada2 */
process q2_dada2 {

    beforeScript '. /appli/bioinfo/qiime/2019.04/env.sh'
    publishDir "${params.outdir}/03_dada2", mode: 'copy', pattern: '*.qz*'
    publishDir "${params.outdir}/00_report", mode: 'copy', pattern: '*_output'

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
    
    //Run only if it is a complete run with full dada2 parameters set
    when :
    !params.first_run

    script :
    """
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
    qiime metadata tabulate \
    --verbose \
    --m-input-file stats.qza \
    --o-visualization stats.qzv >> q2_dada2.log 2>&1
    qiime feature-table summarize \
    --verbose \
    --i-table table.qza \
    --o-visualization table.qzv \
    --m-sample-metadata-file ${metadata} >> q2_dada2.log 2>&1
    qiime feature-table tabulate-seqs \
    --verbose \
    --i-data rep_seqs.qza \
    --o-visualization rep_seqs.qzv >> q2_dada2.log 2>&1
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

    beforeScript '. /appli/bioinfo/qiime/2019.04/env.sh'
    publishDir "${params.outdir}/04_taxonomy", mode: 'copy', pattern: '*.qz*'
    publishDir "${params.outdir}/04_taxonomy", mode: 'copy', pattern: '*.tsv*'
    publishDir "${params.outdir}/00_report", mode: 'copy', pattern: '*_output'

    input :
        file data_repseqs from data_repseqs

    output :
        file 'taxonomy.qza' into data_taxonomy
        file 'taxonomy.qzv' into visu_taxonomy
        file 'ASV_taxonomy.tsv' into taxonomy_tsv
        file 'taxo_output' into taxo_summary

    //Run only if it is a complete run with full dada2 parameters set
    when :
    !params.first_run

    script :
    """
    qiime feature-classifier classify-sklearn \
    --p-n-jobs ${task.cpus} \
    --p-confidence ${params.taxo.confidence} \
    --i-classifier ${params.taxo.database} \
    --i-reads ${data_repseqs} \
    --o-classification taxonomy.qza > q2_taxo.log 2>&1
    qiime metadata tabulate \
    --m-input-file taxonomy.qza \
    --o-visualization taxonomy.qzv >> q2_taxo.log 2>&1
    qiime tools export \
    --input-path taxonomy.qzv \
    --output-path taxo_output >> q2_taxo.log 2>&1
    mv taxo_output/metadata.tsv ASV_taxonomy.tsv >> q2_taxo.log 2>&1
    sed -i '1,2d' ASV_taxonomy.tsv >> q2_taxo.log 2>&1
    sed -i '1 i\\#OTUID\ttaxonomy\tconfidence' ASV_taxonomy.tsv >> q2_taxo.log 2>&1
    """ 
}

/* Prepare output report */
process q2_output {

    beforeScript '. /appli/bioinfo/qiime/2019.04/env.sh'
    publishDir "${params.outdir}/00_report", mode: 'copy', pattern: '*_output'

    input :
        file taxonomy_tsv from taxonomy_tsv
        file dada2_summary from dada2_summary

    output :
        file 'Final_ASV_table_with_taxonomy.biom' into output_biom
        file 'Final_ASV_table_with_taxonomy.tsv' into output_biom_tsv

    //Run only if it is a complete run with full dada2 parameters set
    when :
    !params.first_run

    script :
    """
    biom add-metadata \
    -i ${dada2_summary}/feature-table.biom \
    --observation-metadata-fp ${taxonomy_tsv} \
    -o Final_ASV_table_with_taxonomy.biom \
    --sc-separated taxonomy > q2_output.log 2>&1
    biom convert \
    -i Final_ASV_table_with_taxonomy.biom \
    -o Final_ASV_table_with_taxonomy.tsv \
    --to-tsv \
    --header-key taxonomy >> q2_output.log 2>&1
    """
}

/*
process statisticalanalysis {

    beforeScript ". /appli/bioinfo/R/3.6.1/env.sh" 
    publishDir "${params.outdir}/07_statistical_analysis", mode: 'copy'
    
    input :
        file metadata from metadata
        file Rstatsscript from Rstatsscript
        file output_biom_tsv from output_biom_tsv

    output :
        file 'R/FIGURES' into figures

    //Run only if it is a complete run with full dada2 parameters set
    when :
    !params.first_run
    
    script :
    """
    mkdir -p ${params.stats.data} 
    mkdir -p ${params.stats.script}
    cp ${metadata} ${params.stats.data}
    sed -i 's/#SampleID/SampleID/g' ${params.stats.data}/q2_metadata
    cp ${output_biom_tsv} ${params.stats.data}
    sed -i '1d' ${params.stats.data}/Final_ASV_table_with_taxonomy.tsv
    sed -i 's/#OTU ID/ASV_ID/g' ${params.stats.data}/Final_ASV_table_with_taxonomy.tsv
    Rscript --vanilla ${Rstatsscript} ${projectName} ${params.stats.perc_abund_threshold} ${params.stats.distance} ${params.stats.column_sample_replicat} ${params.stats.exp_var_samples} > stats.log 2&>1
    """
}
*/
