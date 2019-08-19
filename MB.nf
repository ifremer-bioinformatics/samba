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

/* Import metabarcode data */
process q2_import {

    beforeScript '. /appli/bioinfo/qiime/2019.04/env.sh'
    publishDir "${params.outdir}/import", mode: 'copy'

    input : file q2_manifest from manifest

    output : 
        file 'data.qza' into imported_data
        file 'data.qzv' into imported_visu
        file 'summary_output' into imported_summary

    shell :
    """
    qiime tools import \
    --input-path ${q2_manifest} \
    --output-path 'data.qza' \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-format PairedEndFastqManifestPhred33V2 > q2_import.log 2>&1
    qiime demux summarize \
    --verbose \
    --i-data 'data.qza' \
    --o-visualization 'data.qzv' >> q2_import.log 2>&1
    qiime tools export \
    --input-path 'data.qzv' \
    --output-path 'summary_output' >> q2_import.log 2>&1
    """
}

/* Trim metabarcode data with cutadapt */
process q2_cutadapt {

    beforeScript '. /appli/bioinfo/qiime/2019.04/env.sh'
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    
    input : file cutadapt_data from imported_data

    output :
        file 'data_trimmed.qza' into trimmed_data
        file 'data_trimmed.qzv' into trimmed_visu
        file 'trimmed_output' into trimmed_summary

    shell :
    """
    qiime cutadapt trim-paired \
    --verbose \
    --p-cores ${task.cpus} \
    --i-demultiplexed-sequences ${cutadapt_data} \
    --p-front-f ${params.cutadapt.primerF} \
    --p-front-r ${params.cutadapt.primerR} \
    --p-error-rate ${params.cutadapt.errorRate} \
    --p-discard-untrimmed --p-match-read-wildcards \
    --p-overlap ${params.cutadapt.overlap} \
    --o-trimmed-sequences 'data_trimmed.qza' > q2_cutadapt.log 
    qiime demux summarize \
    --verbose \
    --i-data 'data_trimmed.qza' \
    --o-visualization 'data_trimmed.qzv' >> q2_cutadapt.log 2>&1
    qiime tools export \
    --input-path 'data_trimmed.qzv' \
    --output-path 'trimmed_output' >> q2_cutadapt.log 2>&1
    """
}

/* Run dada2 */
process q2_dada2 {

    beforeScript '. /appli/bioinfo/qiime/2019.04/env.sh'
    publishDir "${params.outdir}/dada2", mode: 'copy'

    input : 
        file dada2_data from trimmed_data
        file q2_metadata from metadata        

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
    --i-demultiplexed-seqs ${dada2_data} \
    --p-trim-left-f ${params.dada2.trim3F} \
    --p-trim-left-r ${params.dada2.trim3R} \
    --p-trunc-len-f ${params.dada2.trunclenF} \
    --p-trunc-len-r ${params.dada2.trunclenR} \
    --p-max-ee ${params.dada2.maxee} \
    --p-trunc-q ${params.dada2.minqual} \
    --p-chimera-method ${params.dada2.chimeras} \
    --p-n-threads ${task.cpus} \
    --o-representative-sequences 'rep_seqs.qza' \
    --o-table 'table.qza' \
    --o-denoising-stats 'stats.qza' > q2_dada2.log 2>&1
    qiime metadata tabulate \
    --verbose \
    --m-input-file 'stats.qza' \
    --o-visualization 'stats.qzv' >> q2_dada2.log 2>&1
    qiime feature-table summarize \
    --verbose \
    --i-table 'table.qza' \
    --o-visualization 'table.qzv' \
    --m-sample-metadata-file ${q2_metadata} >> q2_dada2.log 2>&1
    qiime feature-table tabulate-seqs \
    --verbose \
    --i-data 'rep_seqs.qza' \
    --o-visualization 'rep_seqs.qzv' >> q2_dada2.log 2>&1
    qiime tools export \
    --input-path 'rep_seqs.qzv' \
    --output-path 'dada2_output' >> q2_dada2.log 2>&1
    qiime tools export \
    --input-path 'table.qzv' \
    --output-path 'dada2_output' >> q2_dada2.log 2>&1
    qiime tools export \
    --input-path 'stats.qzv' \
    --output-path 'dada2_output' >> q2_dada2.log 2>&1
    """
}
