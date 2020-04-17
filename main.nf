#!/usr/bin/env nextflow

/*
========================================================================================
                         nf-core/samba
========================================================================================
 nf-core/samba Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/samba
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    // Add to this help message with new command line parameters
    log.info nfcoreHeader()
    log.info"""

Usage:

The typical command for running the pipeline is as follows:

	nextflow run main.nf --input_metadata 'PATH-TO-metadata.csv' --input_manifest 'PATH-TO-manifest.csv' --taxo.database "PATH-TO-preformatted-QIIME-db.qza"	-profile conda

	Mandatory arguments:
	--input_metadata		Path to input file with project samples metadata (csv format).
	--input_manifest		Path to input file with samples reads files paths (csv format).
      	-profile			Configuration profile to use. Can use multiple (comma separated).
					Available: conda.
	Generic:
	--singleEnd			Specifies that the input is single-end reads.

	Other options
	--outdir			The output directory where the results will be saved.
	-w/--work-dir			The temporary directory where intermediate data will be saved.
	-name				Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

	Data integrity:
	--data_integrity_enable		Data integrity checking step. Set to false to deactivate this step. (default = true)
	--data_integrity.primerF	Forward primer with '.' characters instead of degenerated bases.
	--data_integrity.primerR	Reverse primer with '.' characters instead of degenerated bases.
	--data_integrity.barcode_filter	Percentage of sample barcode supposed to be found in raw reads (default : 90).
	--data_integrity.primer_filter	Percentage of primers supposed to be found in raw reads (default : 70).
   
	Raw reads cleaning:
	--cutadapt.primerF		Forward primer (to be used in Cutadapt cleaning step).
	--cutadapt.primerR		Reverse primer (to be used in Cutadapt cleaning step).
	--cutadapt.errorRate		Cutadapt error rate allowed to match primers (default : 0.1).
	--cutadapt.overlap		Cutadapt overlaping length between primer and read (default : 18).

	ASVs inference:
	--dada2.trimLeft		The number of nucleotides to remove from the start of each forward read (default : 0 = no trimming).
	--dada2.trimRigth		The number of nucleotides to remove from the start of each reverse read (default : 0 = no trimming).
	--dada2.FtruncLen		Truncate forward reads after FtruncLen bases. Reads shorter than this are discarded (default : 0 = no trimming).
	--dada2.RtruncLen		Truncate reverse reads after RtruncLen bases. Reads shorter than this are discarded (default : 0 = no trimming).
	--dada2.FmaxEE			Forward reads with higher than maxEE "expected errors" will be discarded (default = 2).
	--dada2.RmaxEE			Reverse with higher than maxEE "expected errors" will be discarded (default = 2). 
	--dada2.minQ			After truncation, reads contain a quality score less than minQ will be discarded (default = 10).
	--dada2.chimeras		Chimera detection method : default = "consensus". Set to "pooled" if the samples in the sequence table are all pooled together for bimera identification. 

	Merge ASVs tables:
	--dada2.merge			Set to true to merge Dada2 ASVs tables.
	--dada2.merge_tabledir		Path to the directory containing the ASVs tables to merge (this directory must contain only the ASVs tables to merge).
	--dada2.merge_repseqdir		Path to the directory containing the representative sequences to merge (this directory must constain only the representative sequences to merge).

	Distribution based-clustering:
	--dbotu3_enable			Distribution based-clustering step. Set to false to deactivate this step. (defalut = true)
	--dbotu3.gen_crit		dbotu3 Genetic criterion (default = 0.1).
	--dbotu3.abund_crit		dbotu3 Abundance criterion (default = 10).
	--dbotu3.pval_crit		dbotu3 P-value criterion (default = 0.0005).

	Taxonomic assignation:
	--taxo.extract_db		Set to true to extract specific region from reference database.
	--taxo.seqs_db			Path to reference database.
	--taxo.taxa_db			Path to taxonomic reference database.
	--taxo.database			Path to preformatted QIIME2 format database.
	--taxo.confidence		RDP confidence threshold (default = 90).

	Decontamination:
        --microDecon_enable		Sample decontamination step. Set to true to activate this step. (default = false)
	--microDecon.control_list	Comma separated list of control samples (e.g : "sample1,sample4,sample7").
	--microDecon.nb_controls	Number of controled samples listed.
	--microDecon.nb_samples		Number of samples that are not control samples.
   
	Predict functionnal abundance:
	--picrust2.method		HSP method of your choice. (default = 'mp' ) The most accurate prediction methode. Faster method: 'pic'.
	--picrust2.nsti			Max nsti value accepted. (default = 2) NSTI cut-off of 2 should eliminate junk sequences. 
      
	Statistics:
	--stats_alpha_enable		Set to false to deactivate Alpha diversity statistics step. (default = true)
	--stats_beta_enable		Set to false to deactivate Beta diversity statistics steps. (default = true)
	--stats_sets_analysis_enable	Set to false to deactivate UpsetR graphs steps. (default = true)

	--stats.kingdom			Kingdom to be displayed in barplots.
	--stats.taxa_nb			Number of taxa to be displayed in barplots.
	--stats.alpha_div_group		According to your metadata file, select the column name corresponding to the variable to group samples for Alpha diversity.
	--stats.beta_div_var		According to your metadata file, select the column name corresponding to the variable of interest for Beta diversity.
	--stats.sets_analysis_crit	According to your metadata file, select the column name corresponding to the variable of interest for UpsetR graphs. 
	--stats.hc_method		Hierarchical clustering method (default = 'ward.D2').

	--stats_only			Perform only statistical analysis (ASV table and newick tree required). Set to true to activate. (default = false)
	--inasv_table			if stats_only is activated, set the path to your own ASV table in tsv format.
	--innewick			if stats_only is activated, set the path to your own phylogenetic tree in newick format.
	
	Final analysis report:
	--report_enable			Set to false to deactivate report creation. (default = true)

    """.stripIndent()
}

/**********************************************************
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
if (params.help){
    helpMessage()
    exit 0
}
println "--------------------------------------------------------------"
println "Workflow for project : $params.projectName"
println "Workflow description : $workflow.manifest.description"
println "Workflow gitLab URL : $workflow.manifest.homePage"
println "Workflow authors : $workflow.manifest.author"
println "Workflow source code : $workflow.projectDir"
println "Cmd line: $workflow.commandLine"
println "Workflow working/temp directory : $workflow.workDir"
println "Workflow output/publish directory : $params.outdir"
if (params.singleEnd) println "Going to process SingleEnd data"
println "--------------------------------------------------------------"
println "Optionnal activated steps : "
if (params.stats_only) {
   println "- Statistics steps only activated"
} else {
   if (params.data_integrity_enable) println "- Data integrity step enabled"
   if (params.dbotu3_enable) println "- Distribution based-clustering step enabled"
   if (params.microDecon_enable) println "- Decontamination step enabled"
}
if (params.stats_alpha_enable) println "- Alpha diversity statistics step enabled"
if (params.stats_beta_enable) println "- Beta diversity statistics steps enabled"
if (params.stats_sets_analysis_enable) println "- UpsetR graphs steps enabled"
println "--------------------------------------------------------------"

if(params.dada2.dada2merge == false) {
    Channel.fromPath(params.input_manifest, checkIfExists:true).into { manifest ; manifest4integrity }
    println "Manifest file : $params.input_manifest"
}
Channel.fromPath(params.input_metadata, checkIfExists:true).into { metadata; metadata_dbotu3 ; metadata4stats ; metadata4integrity ; metadata4picrust2 }
println "Metadata file : $params.input_metadata"

//Copy base.config file to output directory for each run
paramsfile = file('conf/base.config')
paramsfile.copyTo("$params.outdir/conf/base.config")

if (params.taxo.extract_db && params.taxo.database == null ) {
   println("ERROR : When extract database option (params.taxo.extract_db) is enable, a taxonomy database (params.taxo.database) must be set.");
   System.exit(1);
}

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

when :
	params.data_integrity_enable && params.stats_only == false && params.dada2.dada2merge == false

script :
"""
${baseDir}/lib/data_integrity.sh ${manifest} ${metadata} ${params.data_integrity.primerF} ${params.data_integrity.primerR} data_integrity.csv verifications.ok verifications.bad ${params.data_integrity.barcode_column_name} ${params.data_integrity.sampleid_column_name} ${params.data_integrity.R1_single_files_column_name} ${params.data_integrity.R1_files_column_name} ${params.data_integrity.R2_files_column_name} ${params.data_integrity.barcode_filter} ${params.data_integrity.primer_filter} ${params.singleEnd} &> data_integrity.log 2>&1
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

	label 'qiime2_env'

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

	when :
		params.stats_only == false && params.dada2.dada2merge == false

	script :
	"""
	${baseDir}/lib/q2_import.sh ${params.singleEnd} ${q2_manifest} data.qza data.qzv import_output completecmd &> q2_import.log 2>&1
	"""
}

/* Trim metabarcode data with cutadapt */
process q2_cutadapt {

	label 'qiime2_env'

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

	when :
		params.stats_only == false && params.dada2.dada2merge == false

	script :
	"""
	${baseDir}/lib/q2_cutadapt.sh ${params.singleEnd} ${task.cpus} ${imported_data} ${params.cutadapt.primerF} ${params.cutadapt.primerR} ${params.cutadapt.errorRate} ${params.cutadapt.overlap} data_trimmed.qza data_trimmed.qzv trimmed_output completecmd &> q2_cutadapt.log 2>&1
	"""
}

/* Run dada2 */
process q2_dada2 {

	label 'qiime2_env' 

	publishDir "${params.outdir}/${params.dada2_dirname}", mode: 'copy', pattern: '*.qz*'
	publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: '*_output'
	publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_dada2 -> "cmd/${task.process}_complete.sh" }

	input : 
		file trimmed_data from trimmed_data
		file metadata from metadata 

	output :
		file 'rep_seqs.qza' into data_repseqs, dada2_seqs_dbotu3, dada2_seqs_taxo, dada2_seqs_decontam, dada2_seqs_phylo, dada2_seqs_picrust2
		file 'rep_seqs.qzv' into visu_repseps
		file 'table.qza' into data_table, dada2_table_picrust2
		file 'table.qzv' into visu_table
		file 'stats.qza' into stats_table
		file 'stats.qzv' into visu_stats
		file 'dada2_output' into dada2_summary
		file 'completecmd' into complete_cmd_dada2

	when :
		params.stats_only == false && params.dada2.dada2merge == false
	
	script :
	"""
	${baseDir}/lib/q2_dada2.sh ${params.singleEnd} ${trimmed_data} ${metadata} rep_seqs.qza rep_seqs.qzv table.qza table.qzv stats.qza stats.qzv dada2_output ${params.dada2.trimLeft} ${params.dada2.trimRigth} ${params.dada2.FtruncLen} ${params.dada2.RtruncLen} ${params.dada2.FmaxEE} ${params.dada2.RmaxEE} ${params.dada2.minQ} ${params.dada2.chimeras} ${task.cpus} completecmd &> q2_dada2.log 2>&1
	"""
}

/* Run dbotu3 */
process q2_dbotu3 {

	label 'qiime2_2019_env'

	publishDir "${params.outdir}/${params.dbotu3_dirname}", mode: 'copy', pattern: '*.qz*'
	publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: '*_output'
	publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_dbotu3 -> "cmd/${task.process}_complete.sh" }

	input :
		file table from data_table
		file seqs from dada2_seqs_dbotu3
		file metadata_dbotu3 from metadata_dbotu3

	output :
		file 'dbotu3_details.txt' into dbotu3_details
		file 'dbotu3_seqs.qza' into dbotu3_seqs_decontam, dbotu3_seqs_taxo, dbotu3_seqs_phylo, dbotu3_seqs_picrust2
		file 'dbotu3_seqs.qzv' into dbotu3_seqs_visu
		file 'dbotu3_table.qza' into dbotu3_table, dbotu3_table_picrust2
		file 'dbotu3_table.qzv' into dbotu3_table_visu
		file 'dbotu3_output' into dbotu3_summary
		file 'completecmd' into complete_cmd_dbotu3

	when :
		params.stats_only == false && params.dada2.dada2merge == false && params.dbotu3_enable

	script :
	"""
	${baseDir}/lib/q2_dbotu3.sh ${table} ${seqs} ${metadata_dbotu3} dbotu3_details.txt dbotu3_seqs.qza dbotu3_seqs.qzv dbotu3_table.qza dbotu3_table.qzv dbotu3_output ${params.dbotu3.gen_crit} ${params.dbotu3.abund_crit} ${params.dbotu3.pval_crit} completecmd &> q2_dbotu3.log 2>&1
	"""
}

/* dada2 merge ASV/seqs */
process q2_dada2_merge {

	label 'qiime2_env'

	publishDir "${params.outdir}/${params.dada2_dirname}/merged", mode: 'copy', pattern: '*.qza'
	publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_dada2merge -> "cmd/${task.process}_complete.sh" }

	input :
		path table_dir from params.dada2.merge_tabledir
		path seq_dir from params.dada2.merge_repseqdir

	output :
		file 'merged_table.qza' into merged_table_picrust2
		file 'merged_seq.qza' into merge_seqs_taxo, merge_seqs_phylo, merge_seqs_picrust2
		file 'merge_output' into merge_summary
		file 'completecmd' into complete_cmd_dada2merge

	when :
	params.dada2.dada2merge && params.stats_only == false

	script :
	"""
	${baseDir}/lib/q2_merge.sh ${table_dir} ${seq_dir} merged_table.qza merged_seq.qza merge_output completecmd &> q2_merge.log 2>&1
	"""
}

if (params.dada2.dada2merge) {
   merge_seqs_taxo.set { seqs_taxo }
   merge_summary.set { summary }
} else {
   if (params.dbotu3_enable) {
      dbotu3_seqs_taxo.set { seqs_taxo }
      dbotu3_summary.into { summary; decontam_summary }
   } else {
      dada2_seqs_taxo.set { seqs_taxo }
      dada2_summary.into { summary; decontam_summary }
   }
}

/* Run taxonomy assignment */
process q2_taxonomy {

	label 'qiime2_env'

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
		file 'ASV_table_with_taxonomy.tsv' into biom_tsv, biom_tsv_decontam
		file 'taxonomic_database.qza' optional true into trained_database
		file 'seqs_db_amplicons.qza' optional true into seqs_db_filtered
		file 'completecmd' into complete_cmd_taxo

	when :
	params.stats_only == false

	script :
	"""
	${baseDir}/lib/q2_taxo.sh ${task.cpus} ${params.taxo.seqs_db} ${params.taxo.taxa_db} ${params.taxo.database} ${params.taxo.extract_db} ${params.cutadapt.primerF} ${params.cutadapt.primerR} ${params.taxo.confidence} ${repseqs_taxo} taxonomy.qza taxonomy.qzv taxo_output ASV_taxonomy.tsv ${summary} ASV_table_with_taxonomy.biom ASV_table_with_taxonomy.tsv taxonomic_database.qza seqs_db_amplicons.qza completecmd &> q2_taxo.log 2>&1
	""" 
}

/* Run sample decontamination using MicroDecon */

process microDecon_step1 {

	label 'microdecon_env'

	publishDir "${params.outdir}/${params.microDecon_dirname}", mode: 'copy', pattern: 'decontaminated_ASV_table.tsv'
	publishDir "${params.outdir}/${params.microDecon_dirname}", mode: 'copy', pattern: 'abundance_removed.txt'
	publishDir "${params.outdir}/${params.microDecon_dirname}", mode: 'copy', pattern: 'ASV_removed.txt'
	publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: 'decontaminated_ASV_table.tsv'
	publishDir "${params.outdir}/${params.report_dirname}/microDecon", mode: 'copy', pattern: 'decontaminated_ASV_table.tsv'
	publishDir "${params.outdir}/${params.report_dirname}/microDecon", mode: 'copy', pattern: 'abundance_removed.txt'
	publishDir "${params.outdir}/${params.report_dirname}/microDecon", mode: 'copy', pattern: 'ASV_removed.txt'
	publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_microDecon -> "cmd/${task.process}_complete.sh" }

	input :
		file microDecon_table from biom_tsv_decontam

	output :
		file 'decontaminated_ASV_table.tsv' into decontam_table, decontam_table_step2, decontam_table_step3
		file 'abundance_removed.txt' into abund_removed
		file 'ASV_removed.txt' into ASV_removed
		file 'completecmd' into complete_cmd_microDecon
		
	when :
		params.stats_only == false && params.dada2.dada2merge == false && params.microDecon_enable == true
		
	shell :
	""" 
	sed '1d' ${microDecon_table} > microDecon_table
	sed -i 's/#OTU ID/ASV_ID/g' microDecon_table
	${baseDir}/lib/microDecon.R microDecon_table ${params.microDecon.control_list} ${params.microDecon.nb_controls} ${params.microDecon.nb_samples} decontaminated_ASV_table.tsv abundance_removed.txt ASV_removed.txt &> microDecon.log 2>&1
	cp ${baseDir}/lib/microDecon.R completecmd &>> microDecon.log 2>&1
	
	"""
}

process microDecon_step2 {

	label 'qiime2_env'

	publishDir "${params.outdir}/${params.microDecon_dirname}", mode: 'copy', pattern: 'decontaminated_ASV_table.qza'
	publishDir "${params.outdir}/${params.report_dirname}/microDecon", mode: 'copy', pattern: 'decontaminated_ASV_table.qza'

	input :
		file table4microDecon from decontam_table_step2

	output :
		file 'decontaminated_ASV_table.qza' into decontam_table_qza, decontam_table_picrust2
		
	when :
		params.stats_only == false && params.dada2.dada2merge == false && params.microDecon_enable == true
		
	shell :
	"""
	biom convert -i ${table4microDecon} -o decontaminated_ASV_table.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy
	qiime tools import --input-path decontaminated_ASV_table.biom --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path decontaminated_ASV_table.qza
	"""
}

/* Extract non-contaminated ASV ID */
process microDecon_step3 {

	label 'seqtk_env'

	publishDir "${params.outdir}/${params.microDecon_dirname}", mode: 'copy', pattern: 'decontaminated_ASV_ID.txt'
	publishDir "${params.outdir}/${params.microDecon_dirname}", mode: 'copy', pattern: 'decontaminated_ASV.fasta'
	publishDir "${params.outdir}/${params.report_dirname}/microDecon", mode: 'copy', pattern: 'decontaminated_ASV.fasta'

	input :
		file decontam_table from decontam_table_step3
		file dada2_summary from decontam_summary

	output :
		file 'decontaminated_ASV_ID.txt' into decontam_ASV_ID
		file 'decontaminated_ASV.fasta' into decontam_ASV_fasta
		
	when :
		params.stats_only == false && params.dada2.dada2merge == false && params.microDecon_enable == true
		
	shell :
	"""
	cut -d \$'\t' -f1 ${decontam_table} | sed '1d' > decontaminated_ASV_ID.txt
	seqtk subseq ${dada2_summary}/sequences.fasta decontaminated_ASV_ID.txt > decontaminated_ASV.fasta
	"""
}

/* Run phylogeny from decontaminated ASV sequences */
process microDecon_step4 {

	label 'qiime2_env'
	
	publishDir "${params.outdir}/${params.phylogeny_dirname}", mode: 'copy', pattern: '*.qza'
	publishDir "${params.outdir}/${params.phylogeny_dirname}", mode: 'copy', pattern: '*.txt'
	publishDir "${params.outdir}/${params.phylogeny_dirname}", mode: 'copy', pattern: 'tree_export_dir'
	publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: 'tree_export_dir'
	publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_phylo -> "cmd/${task.process}_complete.sh" }

	input :
		file ASV_fasta from decontam_ASV_fasta

	output :
		file 'decontam_seqs.qza' into decontam_seqs_qza, decontam_seqs_phylo, decontam_seqs_picrust2
		file 'aligned_repseq.qza' into decontam_aligned_repseq
		file 'masked-aligned_repseq.qza' into decontam_masked_aligned
		file 'tree.qza' into decontam_tree
		file 'tree.log' into decontam_tree_bestmodel_log
		file 'rooted_tree.qza' into decontam_rooted_tree
		file 'tree_export_dir' into decontam_tree_export_dir
		file 'tree_export.log' into decontam_tree_export_log
		file 'tree.nwk' into decontam_newick
		file 'completecmd' into complete_cmd_decontam_phylogeny
		
	when :
		params.stats_only == false && params.dada2.dada2merge == false && params.microDecon_enable == true
		
	shell :
	"""
	qiime tools import --input-path ${ASV_fasta} --output-path decontam_seqs.qza --type 'FeatureData[Sequence]'
	${baseDir}/lib/q2_phylogeny.sh decontam_seqs.qza aligned_repseq.qza masked-aligned_repseq.qza tree.qza tree.log rooted_tree.qza tree_export_dir tree_export.log completecmd &> q2_phylogeny.log 2>&1
	cp tree_export_dir/tree.nwk tree.nwk &>> q2_phylogeny.log 2>&1
	
	"""
}

if (params.dada2.dada2merge) {
   merge_seqs_phylo.set { seqs_phylo }
} else {
   if (params.dbotu3_enable) {
      if (params.microDecon_enable) {
         decontam_seqs_phylo.set { seqs_phylo }
      } else {
         dbotu3_seqs_phylo.set { seqs_phylo }
      }
   } else { 
      dada2_seqs_phylo.set { seqs_pĥylo }
   }
}

/* Run phylogeny construction */
process q2_phylogeny {

	label 'qiime2_env'

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
		file 'tree.nwk' into newick_phylo
		file 'completecmd' into complete_cmd_phylogeny

	when :
	    params.stats_only == false

	script :
	"""
	${baseDir}/lib/q2_phylogeny.sh ${repseqs_phylo} aligned_repseq.qza masked-aligned_repseq.qza tree.qza tree.log rooted_tree.qza tree_export_dir tree_export.log completecmd &> q2_phylogeny.log 2>&1
	cp tree_export_dir/tree.nwk tree.nwk &>> q2_phylogeny.log 2>&1
	"""
}

if (params.dada2.dada2merge){
   merged_table_picrust2.set { table_picrust2 }
   merge_seqs_picrust2.set { seqs_picrust2 }
} else {
   if (params.dbotu3_enable) {
      if (params.microDecon_enable) {
         decontam_table_picrust2.set { table_picrust2 }
         decontam_seqs_picrust2.set { seqs_picrust2 }
      } else {
         dbotu3_table_picrust2.set { table_picrust2 }
         dbotu3_seqs_picrust2.set { seqs_picrust2 }
      }
    } else {
      dada2_table_picrust2.set { table_picrust2 }
      dada2_seqs_picrust2.set { seqs_picrust2 }
   }
}

/* Run functional predictions */
process q2_picrust2_analysis {

	label 'qiime2_2019_env'

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

	when :
		params.stats_only == false

	script :
	"""
	${baseDir}/lib/q2_picrust2.sh ${table_picrust2} ${seqs_picrust2} q2-picrust2_output ${task.cpus} ${params.picrust2.method} ${params.picrust2.nsti} complete_picrust2_cmd &> q2_picrust2.log 2>&1
	"""
}

/* Statistical analysis of functional predictions  */
process q2_picrust2_stats {

    label 'r_stats_env'

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

    when :
        params.stats_only == false

    script :
    """
    ${baseDir}/lib/functional_predictions.R ec_metagenome_predictions_with-descriptions.tsv ko_metagenome_predictions_with-descriptions.tsv pathway_abundance_predictions_with-descriptions.tsv ${metadata4picrust2} ${params.stats.beta_div_var} functional_predictions_NMDS ${params.microDecon_enable} ${params.microDecon.control_list} &> picrust2_stats.log 2>&1
    cp ${baseDir}/lib/functional_predictions.R complete_picrust2_stats_cmd &>> picrust2_stats.log 2>&1
    """
}

if (params.stats_only) {
    Channel.fromPath(params.innewick, checkIfExists:true).set { newick }
    Channel.fromPath(params.inasv_table, checkIfExists:true).set { tsv }
    println "Input ASV table used for statistics steps : $params.inasv_table"
} else {
   if (params.microDecon_enable) decontam_table.set { tsv }
   else biom_tsv.set { tsv }
   newick_phylo.set { newick }
}
	
process prepare_data_for_stats {
    
    label 'r_stats_env'
    
    publishDir "${params.outdir}/${params.report_dirname}/R/DATA", mode: 'copy', pattern : '*.tsv'
    publishDir "${params.outdir}/${params.report_dirname}/R/DATA", mode: 'copy', pattern : '*.rds'
     
    input :
        file metadata from metadata4stats
        file biom_tsv from tsv
        file newick_tree from newick
    
    output :
        file 'ASV_table_with_taxo_for_stats.tsv' into biom_tsv_stats
        file 'metadata_stats.tsv' into metadata_stats, metadata_beta, metadata_beta_rarefied, metadata_beta_deseq2, metadata_beta_css
        file 'phyloseq.rds' into phyloseq_rds, phyloseq_rds_alpha, phyloseq_rds_beta, phyloseq_rds_beta_rarefied, phyloseq_rds_beta_deseq2, phyloseq_rds_beta_css,phyloseq_rds_set
 
    script :
    """
    ${baseDir}/lib/prepare_data_for_stats.sh ${metadata} ${biom_tsv} ASV_table_with_taxo_for_stats.tsv metadata_stats.tsv ${params.microDecon_enable} &> stats_prepare_data.log 2&>1
    Rscript --vanilla ${baseDir}/lib/create_phyloseq_obj.R phyloseq.rds ASV_table_with_taxo_for_stats.tsv metadata_stats.tsv ${params.microDecon_enable} ${params.microDecon.control_list} ${newick_tree} &>> stats_prepare_data.log 2&>1 
    """
}
   
process stats_alpha {
   
    label 'r_stats_env'
   
    publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/alpha_diversity", mode: 'copy', pattern : 'index_significance_tests.txt'
    publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/alpha_diversity/diversity_index", mode: 'copy', pattern : 'alpha_div_plots*'
    publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/alpha_diversity/diversity_barplots", mode: 'copy', pattern : 'barplot_*'
    publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/alpha_diversity", mode: 'copy', pattern : 'rarefaction_curve*'
    
    input :
        file phyloseq_rds from phyloseq_rds_alpha

    output :
        file 'alpha_div_plots*' into alpha_div_plots
        file 'index_significance_tests.txt' into index_significance_tests
        file 'barplot_phylum*'into barplot_phylum
        file 'barplot_class*' into barplot_class
        file 'barplot_order*' into barplot_order
        file 'barplot_family*'into barplot_family
        file 'barplot_genus*' into barplot_genus
        file 'rarefaction_curve*' into rarefaction_curve

    when :
        params.stats_alpha_enable
    
    shell :
    """
    Rscript --vanilla ${baseDir}/lib/alpha_diversity.R phyloseq.rds ${params.stats.distance} alpha_div_plots ${params.stats.kingdom} ${params.stats.taxa_nb} barplot_phylum barplot_class barplot_order barplot_family barplot_genus ${params.stats.alpha_div_group} index_significance_tests.txt $workflow.projectDir rarefaction_curve &> stats_alpha_diversity.log 2>&1
    """
}

process stats_beta {

    label 'r_stats_env'

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

    when :
        params.stats_beta_enable

    shell :
    """
    Rscript --vanilla ${baseDir}/lib/beta_diversity.R ${phyloseq_rds} ${params.stats.beta_div_var} ${metadata} $workflow.projectDir NMDS_ PCoA_ ${params.stats.hc_method} hclustering_ variance_significance_tests_ pie_ExpVar_ &> stats_beta_diversity.log 2>&1
    """
}

process stats_beta_rarefied {

    label 'r_stats_env'

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

    when :
        params.stats_beta_enable

    shell :
    """
    Rscript --vanilla ${baseDir}/lib/beta_diversity_rarefied.R ${phyloseq_rds} Final_rarefied_ASV_table_with_taxonomy.tsv ${params.stats.beta_div_var} ${metadata} $workflow.projectDir NMDS_rarefied_ PCoA_rarefied_ ${params.stats.hc_method} hclustering_rarefied_ variance_significance_tests_rarefied_ pie_ExpVar_rarefied_ &> stats_beta_diversity_rarefied.log 2>&1
    """
}

process stats_beta_deseq2 {

    label 'r_stats_env'

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
        
    when :
        params.stats_beta_enable

    shell :
    """
    Rscript --vanilla ${baseDir}/lib/beta_diversity_deseq2.R ${phyloseq_rds} Final_DESeq2_ASV_table_with_taxonomy.tsv ${params.stats.beta_div_var} ${metadata} $workflow.projectDir NMDS_DESeq2_ PCoA_DESeq2_ ${params.stats.hc_method} hclustering_DESeq2_ variance_significance_tests_DESeq2_ pie_ExpVar_DESeq2_ &> stats_beta_diversity_deseq2.log 2>&1
    """
}

process stats_beta_css {

    label 'r_stats_env'

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


    when :
        params.stats_beta_enable

    shell :
    """
    Rscript --vanilla ${baseDir}/lib/beta_diversity_css.R ${phyloseq_rds} Final_CSS_ASV_table_with_taxonomy.tsv ${params.stats.beta_div_var} ${metadata} $workflow.projectDir NMDS_CSS_ PCoA_CSS_ ${params.stats.hc_method} hclustering_CSS_ variance_significance_tests_CSS_ pie_ExpVar_CSS_ &> stats_beta_diversity_css.log 2>&1
    """
}

process stats_sets_analysis {

    label 'r_stats_env'

    publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/", mode: 'copy', pattern : 'upset_plot*'

    input :
        file phyloseq_rds from phyloseq_rds_set

    output :
        file 'upset_plot*' into upset_plot
        file 'end_analysis.ok' into end_ok

    when :
        params.stats_sets_analysis_enable

    shell :
    """
    Rscript --vanilla ${baseDir}/lib/sets_analysis.R ${phyloseq_rds} ${params.stats.sets_analysis_crit} upset_plot &> stats_sets_analysis.log 2>&1
    touch end_analysis.ok
    """
}

if (params.report_enable) {
    Channel.fromPath(params.reportHTML, checkIfExists:true).set { reportHTML }
    Channel.fromPath(params.reportMD, checkIfExists:true).set { reportMD }
}
process report {

    publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'Report_*'

    input :
        file reportHTML from reportHTML
        file reportMD from reportMD
        file 'end_analysis.ok' from end_ok

    output :
        file 'Report_*' into Reports

    when :
       params.report_enable

    shell :
    """
    cp ${reportHTML} Report_${params.projectName}.html
    cp ${reportMD} Report_${params.projectName}.md
    """
}


/* Other functions */
def nfcoreHeader() {
    // Log colors ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/samba v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

