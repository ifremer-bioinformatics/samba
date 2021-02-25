#!/usr/bin/env nextflow
/*
========================================================================================
                         samba
========================================================================================
 samba Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/ifremer-bioinformatics/samba
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    // Add to this help message with new command line parameters
    log.info SeBiMERHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline after filling the conf/custom.config file is as follows:

	nextflow run main.nf -profile <docker,singularity,conda>,custom

	Mandatory arguments:
	--input_metadata [file]		Path to input file with project samples metadata (csv format).
	--input_manifest [file]		Path to input file with samples reads files paths (csv format).
	-profile [str]			Configuration profile to use. Can use multiple (comma separated).
					Available: conda, docker, singularity, test, custom.
	Generic:
	--singleEnd [bool]		Set to true to specify that the inputs are single-end reads.
	--longreads [bool]		Set to true to specify that the inputs are long reads (Nanopore/Pacbio) (default = false for illumina short reads).
    --compress_result [bool]        Zip the final result directory (default = true) 

	Other options
	--outdir [path]			The output directory where the results will be saved.
	-w/--work-dir			The temporary directory where intermediate data will be saved.
	--email [email]                 Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
	--email_on_fail [email]         Same as --email, except only send mail if the workflow is not successful
	-name [str]			Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
	--projectName [str]		Name of the project being analyzed.

	Data integrity:
	--data_integrity_enable [bool]	Data integrity checking step. Set to false to deactivate this step. (default = true)
	--primer_filter [str]		Percentage of primers supposed to be found in raw reads (default : 70).

	Raw reads cleaning:
	--primerF [str]			Forward primer (to be used in Cutadapt cleaning step).
	--primerR [str]			Reverse primer (to be used in Cutadapt cleaning step).
	--errorRate [str]		Cutadapt error rate allowed to match primers (default : 0.1).
	--overlap [str]			Cutadapt overlaping length between primer and read (default : 18 for test dataset, must be changed for user dataset).

	ASVs inference:
	--FtrimLeft [str]		The number of nucleotides to remove from the start of each forward read (default : 0 = no trimming).
	--RtrimLeft [str]		The number of nucleotides to remove from the start of each reverse read (default : 0 = no trimming).
	--FtruncLen [str]		Truncate forward reads after FtruncLen bases. Reads shorter than this are discarded (default : 0 = no trimming).
	--RtruncLen [str]		Truncate reverse reads after RtruncLen bases. Reads shorter than this are discarded (default : 0 = no trimming).
	--FmaxEE [str]			Forward reads with higher than maxEE "expected errors" will be discarded (default = 2).
	--RmaxEE [str]			Reverse with higher than maxEE "expected errors" will be discarded (default = 2).
	--minQ [str]			Truncate reads at the first instance of a quality score less than or equal to minQ (default = 2).
	--chimeras [str]		Chimera detection method : default = "consensus". Set to "pooled" if the samples in the sequence table are all pooled together for bimera identification.

	Merge ASVs tables:
	--dada2merge [bool]		Set to true to merge DADA2 ASVs tables.
	--merge_tabledir [path]		Path to the directory containing the ASVs tables to merge (this directory must contain only the ASVs tables to merge).
	--merge_repseqsdir [path]	Path to the directory containing the representative sequences to merge (this directory must constain only the representative sequences to merge).

	Distribution based-clustering:
	--dbotu3_enable	[bool]		Distribution based-clustering step. Set to false to deactivate this step. (default = true)
	--gen_crit [str]		dbotu3 Genetic criterion (default = 0.1).
	--abund_crit [str]		dbotu3 Abundance criterion (default = 10).
	--pval_crit [str]		dbotu3 P-value criterion (default = 0.0005).

	Taxonomic assignation:
	--extract_db [bool]		Set to true to extract specific region from reference database (default = false)
	--seqs_db [file]		Path to reference database (required if extract_db = true).
	--taxo_db [file]		Path to taxonomic reference database (required if extract_db = true).
	--database [file]		Path to preformatted QIIME2 format database (required if extract_db = false).
	--confidence [str]		Confidence threshold for limiting taxonomic depth. Set to "disable" to disable confidence calculation, or 0 to calculate confidence but not apply it to limit the taxonomic depth of the assignments (default = 0.9).

	Taxonomy filtering:
	--filtering_tax_enable [bool]	Set to true to filter asv table and sequences based on taxonomic assignation (default = false)
	--tax_to_exclude [str]		List of taxa you want to exclude (comma-separated list).

	Decontamination:
	--microDecon_enable [bool]	Sample decontamination step. Set to true to activate this step. (default = false)
	--control_list [str]		Comma separated list of control samples (e.g : "sample1,sample4,sample7") (required if microDecon_enable = true).
	--nb_controls [str]		Number of control sample listed (required if microDecon_enable = true).
	--nb_samples [str]		Number of samples that are not control samples (required if microDecon_enable = true).

	Predict functionnal abundance:
	--picrust2_enable [bool]	Set to true to enable functionnal prediction step. (default = false)
        --picrust_var [str]		According to your metadata file, list the column names corresponding to the variables to group samples for functional predictions (comma-separated list).
	--method [str]			HSP method of your choice. (default = 'mp' ) The most accurate prediction methode. Faster method: 'pic'.
	--nsti [str]			Max NSTI value accepted. (default = 2) NSTI cut-off of 2 should eliminate junk sequences.

	Differential abundance testing:
	--ancom_var [str]	        According to your metadata file, list the column names corresponding to the variables to group samples for ANCOM analysis (comma-separated list).

	Statistics:
	--stats_alpha_enable [bool]	Set to false to deactivate Alpha diversity statistics step. (default = true)
	--stats_beta_enable [bool]	Set to false to deactivate Beta diversity statistics steps. (default = true)
	--stats_desc_comp_enable [bool]	Set to false to deactivate Descriptive comparisons steps. (default = true)

	--kingdom [str]			Kingdom to be displayed in barplots.
	--taxa_nb [str]			Number of taxa to be displayed in barplots.
	--alpha_div_group [str]	 	According to your metadata file, list the column names corresponding to the variables to group samples for Alpha diversity (comma-separated list).
	--beta_div_var [str]		According to your metadata file, select the column name corresponding to the variable of interest for Beta diversity (comma-separated list).
	--desc_comp_crit [str]		According to your metadata file, select the column name corresponding to the variable of interest for Descriptive comparisons graphs (comma-separated list).
	--hc_method [str]		Hierarchical clustering method (default = 'ward.D2').

	--stats_only [bool]		Perform only statistical analysis (ASV table and newick tree required). Set to true to activate. (default = false)
	--inasv_table [file]		if stats_only is activated, set the path to your own ASV table in tsv format.
	--innewick [file]		if stats_only is activated, set the path to your own phylogenetic tree in newick format.

	Final analysis report:
	--report_enable	[bool]		Set to false to deactivate report creation. (default = true)

	Long reads options:
	--lr_type [str]			Long reads technology. For pacbio, [map-pb] and for nanopore, [map-ont]
	--lr_tax_fna [file]		Path to reference database indexed with Minimap2 (required).
	--lr_taxo_flat [file]		Path to taxonomic reference file (required).
	--lr_rank [int]			Minimal rank level to keep a hit as assigned [5]. 1:Kingdom, 2:Phylum, 3:Class, 4:Order, 5:Family, 6:Genus, 7:Species

    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*
 * SET UP CONFIGURATION VARIABLES
 */
import java.text.SimpleDateFormat
def d = new Date()
def dtFormat = "dd/MM/yyyy"
def date = d.format(dtFormat, TimeZone.getTimeZone('GMT+2'))

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

//Copy config files to output directory for each run
paramsfile = file("$baseDir/conf/base.config", checkIfExists: true)
paramsfile.copyTo("$params.outdir/conf/base.config")

if (workflow.profile.contains('shortreadstest')) {
   testparamsfile = file("$baseDir/conf/shortreadstest.config", checkIfExists: true)
   testparamsfile.copyTo("$params.outdir/conf/shortreadstest.config")
}
if (workflow.profile.contains('longreadstest')) {
   testparamsfile = file("$baseDir/conf/longreadstest.config", checkIfExists: true)
   testparamsfile.copyTo("$params.outdir/conf/longreadstest.config")
}
if (workflow.profile.contains('custom')) {
   customparamsfile = file("$baseDir/conf/custom.config", checkIfExists: true)
   customparamsfile.copyTo("$params.outdir/conf/custom.config")
}

//If only running the stats processes of the pipeline
if (params.stats_only) {
   if (!params.inasv_table || params.inasv_table.isEmpty()) {
      log.error "Parameter --inasv_table cannot be null or empty. Set the path to the ASV count table"
      exit 1
   } else {
      inasv_table_ch = Channel.fromPath(params.inasv_table, checkIfExists:true)
                           .set { tsv_only }
   }
   if (!params.innewick || params.innewick.isEmpty()) {
      log.error "Parameter --innewick cannot be null or empty. Set the path to the newick tree"
      exit 1
   } else {
      newick_ch = Channel.fromPath(params.innewick, checkIfExists:true)
                         .set { newick_only }
   }

   //Force to false other processes options
   params.dada2merge = false
   params.data_integrity_enable = false
   params.dbotu3_enable = false
   params.filtering_tax_enable = false
   params.microDecon_enable = false
   params.picrust2_enable = false
} else {
   params.inasv_table = ""
   params.innewick = ""
   inasv_table_ch = Channel.empty()
   newick_ch = Channel.empty()
}

//If pipeline runs for merging dada2 data
if (params.dada2merge) {
   if (!params.merge_tabledir || params.merge_tabledir.isEmpty()) {
      log.error "Parameter --merge_tabledir cannot be null or empty. Set the path to the folder with ASV tables to merge"
      exit 1
   } else {
      dada2merge_tabledir_ch = Channel.fromPath(params.merge_tabledir, checkIfExists:true)
   }
   if (!params.merge_repseqsdir || params.merge_repseqsdir.isEmpty()) {
      log.error "Parameter --merge_repseqsdir cannot be null or empty. Set the path to the folder with ASV sequences to merge"
      exit 1
   } else {
      dada2merge_repseqsdir_ch = Channel.fromPath(params.merge_repseqsdir, checkIfExists:true)
   }
   //Force to false other processes options
   params.data_integrity_enable = false
   params.dbotu3_enable = false
   params.filtering_tax_enable = false
   params.microDecon_enable = false
   params.picrust2_enable = false
   params.stats_only = false
} else {
   params.merge_tabledir = ""
   params.merge_repseqsdir = ""
   dada2merge_tabledir_ch = Channel.empty()
   dada2merge_repseqsdir_ch = Channel.empty()
}

//If taxonomy step require to extract primer regions
if (!params.longreads && !params.stats_only) {
   if (params.extract_db) {
      if (!params.seqs_db || params.seqs_db.isEmpty()) {
         log.error "Parameter --seqs_db cannot be null or empty. Set the path to the reference sequences"
         exit 1
      } 
      if (!params.taxo_db || params.taxo_db.isEmpty()) {
         log.error "Parameter --taxo_db cannot be null or empty. Set the path to the reference taxonomy"
         exit 1
      }
   } else {
      if (!params.database || params.database.isEmpty()) {
         log.error "Parameter --database cannot be null or empty. Set the path to the reference database"
      }
   }
}

if (params.microDecon_enable && !params.control_list) {
   log.error "ERROR : A comma separated list of control samples (--control_list) must be set."
   exit 1
}

if (params.longreads) {
//Force to false other processes options
   params.data_integrity_enable = false
   params.qiime2 = false
   params.dada2merge = false
   params.dbotu3_enable = false
   params.filtering_tax_enable = false
   params.microDecon_enable = false
   params.ancom_enable = false
   params.asv = false
   params.picrust2_enable = false
// Initialized variables not used in long reads analysis
   params.primerF = ""
   params.primerR = ""
   params.errorRate = ""
   params.overlap = ""
   params.FtrimLeft = ""
   params.RtrimLeft = ""
   params.FtruncLen = ""
   params.RtruncLen = ""
   params.FmaxEE = ""
   params.RmaxEE = ""
   params.minQ = ""
   params.chimeras = ""
   params.database = ""
   params.seqs_db = ""
   params.taxo_db = ""
   params.tax_to_exclude = ""
   params.nb_controls = ""
   params.nb_samples = ""
   params.remove_sample = false
   params.sample_to_remove = ""
   params.ancom_var = ""
} else {
// Initialized variables not used in short reads analysis
   params.lr_rank = ""
   params.lr_type = ""
   params.lr_tax_fna = ""
   params.lr_taxo_flat = ""
}

if (!params.stats_alpha_enable) {
   params.alpha_div_group = ""
   params.kingdom = ""
   params.taxa_nb = ""
}

if (!params.stats_beta_enable) {
   params.beta_div_var = ""
   params.hc_method = ""
}

if (!params.desc_comp_enable) {
   params.desc_comp_crit = ""
   params.desc_comp_tax_level = ""
}

if (!params.filtering_tax_enable) {
   params.tax_to_exclude = ""
}

if (!params.remove_sample) {
   params.sample_to_remove = ""
}

if (!params.extract_db) {
   params.seqs_db = ""
   params.taxo_db = ""
}

/*
 * PIPELINE INFO
 */

// Header log info
log.info SeBiMERHeader()
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
summary['Project Name']     = params.projectName
summary['Manifest File']    = params.input_manifest
summary['Metadata File']    = params.input_metadata
summary['Data Type']        = params.singleEnd ? 'Single-End' : 'Paired-End'
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
summary['Config Profile'] = workflow.profile

if (params.stats_only) summary['Stats only'] = "Pipeline running only statistics processes"
if (params.data_integrity_enable) summary['Data integrity'] = "Data integrity checking process enabled"
if (params.dbotu3_enable) summary['Clustering'] = "Distribution based-clustering process enabled"
if (params.filtering_tax_enable) summary['Tax filtering'] = "Filtering ASV table and sequences based on taxonomic assignation"
if (params.microDecon_enable) summary['Decontamination'] = "Sample decontamination process enabled"
if (params.dada2merge) summary['Dada2 merge'] = "Dada2 merge process enabled"
if (params.picrust2_enable) summary['Prediction'] = "Functionnal prediction process enabled"
if (params.ancom_enable) summary['Ancom'] = "Differential abundance ANCOM analysis process enabled"
if (params.stats_alpha_enable) summary['Alpha div'] = "Alpha diversity indexes process enabled"
if (params.stats_beta_enable) summary['Beta div'] = "Beta diversity statistics processes enabled"
if (params.stats_desc_comp_enable) summary["Desc comp"] = "Descriptive comparisons statistics process enabled"
if (params.longreads) summary["Long reads"] = "Going to run long reads processes"

if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

// Check the hostnames against configured profiles
checkHostname()

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'samba-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'samba Workflow Summary'
    section_href: 'https://github.com/ifremer-bioinformatics/samba'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }

// Check the hostnames against configured profiles
checkHostname()

/*
 * STEP 0 - Get test data if running in test profile
*/
if (workflow.profile.contains('test')) {
   process get_test_data {
      label 'internet_access'
      output :
         file 'data_is_ready' into ready_integrity, ready_import, ready_lr
         file 'manifest' into testmanifest
         file 'metadata' into testmetadata
       when :
         !params.stats_only && !params.dada2merge
      script :
      def datatype = params.longreads ? "longreads" : "shortreads"
      """
      get_test_data.sh ${baseDir} 'data_is_ready' ${datatype} 'manifest' 'metadata' &> get_test_data.log 2>&1
      """
   }

   if (!params.stats_only || !params.dada2merge) {
      if (params.longreads) {
        testmetadata.set { metadata4stats } 
      } else {
        testmanifest.into { manifest ; manifest4integrity }
        testmetadata.into { metadata4dada2 ; metadata4dbotu3 ; metadata_filtering_tax ; metadata4stats ; metadata4integrity ; metadata4picrust2 ; metadata4ancom }
      } 
   }
  
} else {
   data_ready = Channel.value("none")
   data_ready.into { ready_integrity ; ready_import ; ready_lr}

   if (!params.stats_only || !params.dada2merge) {
      // Check if input metadata file exits
      if (!params.input_metadata || params.input_metadata.isEmpty()) {
         log.error "Parameter --input_metadata cannot be null or empty. Set the path to the Metadata file."
         exit 1
      } else {
         Channel.fromPath(params.input_metadata, checkIfExists:true)
                .into { metadata4dada2 ; metadata4dbotu3 ; metadata_filtering_tax ; metadata4stats ; metadata4integrity ; metadata4picrust2 ; metadata4ancom }
      }
      if (!params.input_manifest || params.input_manifest.isEmpty()) {
         log.error "Parameter --input_manifest cannot be null or empty. Set the path to the Manifest file."
         exit 1
      } else {
         Channel.fromPath(params.input_manifest, checkIfExists:true)
                .into { manifest ; manifest4integrity }
      }
   } else {
      Channel.empty().into { manifest ; manifest4integrity }
   }
}

/*
 * STEP 1 - Checking data integrity
 */
if (!params.longreads) {
    if (params.data_integrity_enable) {
        /* Check data integrity */
        process data_integrity {
    
            label 'biopython_env'
    
        	publishDir "${params.outdir}/${params.data_integrity_dirname}", mode: 'copy', pattern: 'data_integrity.txt'
        	publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: 'data_integrity.txt'
        	publishDir "${params.outdir}/${params.data_integrity_dirname}", mode: 'copy', pattern: '*.sort'
        	publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: '*.sort'
    
        	input :
        		file manifest from manifest4integrity
        		file metadata from metadata4integrity
                        file ready from ready_integrity
    
        	output :
        		file 'data_integrity.txt'
        		file "${metadata}.sort" into metadata_sort
        		file "${manifest}.sort" into manifest_sort
    
        	when :
        		params.data_integrity_enable && !params.stats_only && !params.dada2merge
        	script :
        	def datatype = params.singleEnd ? "single" : "paired"
        	"""
                data_integrity.py -e ${metadata} -a ${manifest} -p ${params.primer_filter} -r ${datatype} -t ${task.cpus} -c ${params.control_list} &> data_integrity.log 2>&1
                """
        }
        metadata_sort.into { metadata4dada2 ; metadata4dbotu3 ; metadata4stats ; metadata4picrust2 ; metadata4ancom }
        manifest_sort.set { manifest }
    }
    
    /*
     * STEP 2 - Import metabarcoding data into QIIME2 object
     */
    process q2_import {
    
    	label 'qiime2_env1cpu'
    
    	publishDir "${params.outdir}/${params.import_dirname}", mode: 'copy', pattern: 'data.qz*'
    	publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: '*_output'
    	publishDir "${params.outdir}/${params.report_dirname}/version", mode: 'copy', pattern: 'v_*.txt'
    	publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_import -> "cmd/${task.process}_complete.sh" }
    
    	input :
    		file q2_manifest from manifest
                file ready from ready_import
    
    	output :
    		file 'data.qza' into imported_data
    		file 'data.qzv' into imported_visu
    		file 'import_output' into imported_output
    		file 'completecmd' into complete_cmd_import
                    file 'v_qiime2.txt' into qiime2_version
    
    	when :
    		!params.stats_only && !params.dada2merge
    
    	script :
    	"""
    	q2_import.sh ${params.singleEnd} ${q2_manifest} data.qza data.qzv import_output completecmd &> q2_import.log 2>&1
    	qiime --version|grep 'q2cli'|cut -d' ' -f3 > v_qiime2.txt
    	"""
    }
    
    
    /*
     * STEP 3 - Trim metabarcode data with cutadapt
     */
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
    		file 'trimmed_output' into trimmed_output
    		file 'completecmd' into complete_cmd_cutadapt
    
    	when :
    		!params.stats_only && !params.dada2merge
    
    	script :
    	"""
    	q2_cutadapt.sh ${params.singleEnd} ${task.cpus} ${imported_data} ${params.primerF} ${params.primerR} ${params.errorRate} ${params.overlap} data_trimmed.qza data_trimmed.qzv trimmed_output completecmd &> q2_cutadapt.log 2>&1
    	"""
    }
    
    /*
     * STEP 4 - Dada2 ASVs inference
     */
    process q2_dada2 {
    
    	label 'qiime2_env'
    
    	publishDir "${params.outdir}/${params.dada2_dirname}", mode: 'copy', pattern: '*.qz*'
    	publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: '*_output'
    	publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_dada2 -> "cmd/${task.process}_complete.sh" }
    
    	input :
    		file trimmed_data from trimmed_data
    		file metadata from metadata4dada2
    
    	output :
    		file 'rep_seqs.qza' into dada2_seqs_dbotu3, dada2_seqs_taxo, dada2_seqs_filtering_tax, dada2_seqs_decontam, dada2_seqs_phylo, dada2_seqs_picrust2, dada2_seqs_ancom
    		file 'rep_seqs.qzv' into visu_repseps
    		file 'table.qza' into dada2_table_dbotu3, dada2_table_filtering_tax, dada2_table_picrust2, dada2_table_ancom
    		file 'table.qzv' into visu_table
    		file 'stats.qza' into stats_table
    		file 'stats.qzv' into visu_stats
    		file 'dada2_output' into dada2_output
    		file 'completecmd' into complete_cmd_dada2
    
    	when :
    		!params.stats_only && !params.dada2merge
    
    	script :
    	"""
    	q2_dada2.sh ${params.singleEnd} ${trimmed_data} ${metadata} rep_seqs.qza rep_seqs.qzv table.qza table.qzv stats.qza stats.qzv dada2_output ${params.FtrimLeft} ${params.RtrimLeft} ${params.FtruncLen} ${params.RtruncLen} ${params.FmaxEE} ${params.RmaxEE} ${params.minQ} ${params.chimeras} ${task.cpus} completecmd &> q2_dada2.log 2>&1
    	"""
    }
    
    /*
     * STEP 5 - Distribution based-clustering with dbotu3
     */
    if (params.dbotu3_enable) {
        /* Run dbotu3 */
        process q2_dbotu3 {
    
          label 'qiime2_env'
    
        	publishDir "${params.outdir}/${params.dbotu3_dirname}", mode: 'copy', pattern: '*.qz*'
        	publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: '*_output'
        	publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_dbotu3 -> "cmd/${task.process}_complete.sh" }
    
        	input :
        		file table from dada2_table_dbotu3
        		file seqs from dada2_seqs_dbotu3
        		file metadata from metadata4dbotu3
    
        	output :
        		file 'dbotu3_details.txt' into dbotu3_details
        		file 'dbotu3_seqs.qza' into dbotu3_seqs_decontam, dbotu3_seqs_taxo, dbotu3_seqs_filtering_tax, dbotu3_seqs_phylo, dbotu3_seqs_picrust2, dbotu3_seqs_ancom
        		file 'dbotu3_seqs.qzv' into dbotu3_seqs_visu
        		file 'dbotu3_table.qza' into dbotu3_table, dbotu3_table_filtering_tax, dbotu3_table_picrust2, dbotu3_table_ancom
        		file 'dbotu3_table.qzv' into dbotu3_table_visu
        		file 'dbotu3_output' into dbotu3_output
        		file 'completecmd' into complete_cmd_dbotu3
    
        	when :
        		!params.stats_only && !params.dada2merge && params.dbotu3_enable
    
        	script :
        	"""
        	q2_dbotu3.sh ${table} ${seqs} ${metadata} dbotu3_details.txt dbotu3_seqs.qza dbotu3_seqs.qzv dbotu3_table.qza dbotu3_table.qzv dbotu3_output ${params.gen_crit} ${params.abund_crit} ${params.pval_crit} completecmd &> q2_dbotu3.log 2>&1
        	"""
        }
    }
    
    /*
     * STEP 6 - Use Dada2 merge to merge ASVs tables and sequences
     */
    if (params.dada2merge) {
        process q2_dada2_merge {
    
                label 'qiime2_env'
    
                publishDir "${params.outdir}/${params.dada2_dirname}/merged", mode: 'copy', pattern: '*.qza'
                publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_dada2merge -> "cmd/${task.process}_complete.sh" }
    
                input :
                        path table_dir from dada2merge_tabledir_ch
                        path seq_dir from dada2merge_repseqsdir_ch
    
                output :
                        file 'merged_table.qza' into merge_table_picrust2, merge_table_ancom
                        file 'merged_seq.qza' into merge_seqs_taxo, merge_seqs_phylo, merge_seqs_picrust2, merge_seqs_ancom
                        file 'merge_output' into merge_output
                        file 'completecmd' into complete_cmd_dada2merge
    
                when :
                    params.dada2merge
    
                script :
                """
                q2_merge.sh ${table_dir} ${seq_dir} merged_table.qza merged_seq.qza merge_output completecmd &> q2_merge.log 2>&1
                """
        }
    }
    
    outputA = params.dada2merge ? merge_output : dada2_output
    output_ch = params.dbotu3_enable ? dbotu3_output : outputA
    output_ch.into { taxonomy_output ; decontam_output }
    
    seqs_taxoA = params.dada2merge ? merge_seqs_taxo : dada2_seqs_taxo
    seqs_taxo = params.dbotu3_enable ? dbotu3_seqs_taxo : seqs_taxoA
   
    /*
    * STEP 7  -  Run taxonomy assignment
    */
    process q2_taxonomy {
           label 'qiime2_highRAM'

           publishDir "${params.outdir}/${params.taxo_dirname}", mode: 'copy', pattern: '*.qz*'
           publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: 'taxo_output'
           publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: 'ASV_taxonomy.tsv'
           publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: 'ASV_table*'
           publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_taxo -> "cmd/${task.process}_complete.sh" }

           input : 
                file repseqs_taxo from seqs_taxo
    		file summary_output from taxonomy_output
           output :
    		file 'taxonomy.qza' into data_taxonomy_filtering_tax, data_taxonomy_ancom
    		file 'taxonomy.qzv' into visu_taxonomy
    		file 'ASV_taxonomy.tsv' into taxonomy_tsv, taxonomy_tsv_filtering_tax
    		file 'taxo_output' into taxo_output
    		file 'ASV_table_with_taxonomy.biom' into biom
    		file 'ASV_table_with_taxonomy.tsv' into biom_tsv, biom_tsv_decontam
    		file 'taxonomic_database.qza' optional true into trained_database
    		file 'seqs_db_amplicons.qza' optional true into seqs_db_filtered
    		file 'completecmd' into complete_cmd_taxo
           when :
                !params.stats_only

           script :
           """
           q2_taxo.sh ${task.cpus} ${params.extract_db} ${params.primerF} ${params.primerR} ${params.confidence} ${repseqs_taxo} taxonomy.qza taxonomy.qzv taxo_output ASV_taxonomy.tsv ${summary_output} ASV_table_with_taxonomy.biom ASV_table_with_taxonomy.tsv taxonomic_database.qza seqs_db_amplicons.qza completecmd tmpdir ${params.database} ${params.seqs_db} ${params.taxo_db} &> q2_taxo.log 2>&1
           """ 
    }

    asv_table = params.dbotu3_enable ? dbotu3_table_filtering_tax : dada2_table_filtering_tax
    asv_seq = params.dbotu3_enable ? dbotu3_seqs_filtering_tax : dada2_seqs_filtering_tax

     /*
     * STEP 8 -  Filtering ASVs based on taxonomy assignment
     */
     if (params.filtering_tax_enable) {
         process q2_filtering_tax {
        
             label 'qiime2_env'

             publishDir "${params.outdir}/${params.tax_filtering_dirname}", mode: 'copy', pattern: '*.qz*'
             publishDir "${params.outdir}/${params.report_dirname}/tax_filtering", mode: 'copy', pattern: 'tax_filtered_table_with_tax.*'
             publishDir "${params.outdir}/${params.report_dirname}/tax_filtering", mode: 'copy', pattern: 'tax_filtered_output'
             publishDir "${params.outdir}/${params.report_dirname}/version", mode: 'copy', pattern: 'v_*.txt'
             publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_filtering_tax -> "cmd/${task.process}_complete.sh" }

             input :
                 file asv_table from asv_table
                 file asv_tax from data_taxonomy_filtering_tax
                 file asv_seq from asv_seq
                 file asv_tax_tsv from taxonomy_tsv_filtering_tax
                 file metadata from metadata_filtering_tax
             output :
                 file 'tax_filtered_table.qza' into tax_filtered_table, tax_filtered_table_picrust2, tax_filtered_table_ancom
                 file 'tax_filtered_seq.qza' into tax_filtered_seq, tax_filtered_seq_picrust2
                 file 'tax_filtered_table.qzv' into tax_filtered_table_visu
                 file 'tax_filtered_output' into tax_filtered_output
                 file 'tax_filtered_seq.qzv' into tax_filtered_seq_visu
                 file 'tax_filtered_table_with_tax.biom' into tax_filtered_table_biom
                 file 'tax_filtered_table_with_tax.tsv' into tax_filtered_table_tsv, tax_filtered_table_tsv_stats
                 file 'completecmd' into complete_cmd_filtering_tax
             when :
                 !params.stats_only
             script :
             """
             q2_filtering_tax.sh ${asv_table} ${asv_tax} ${params.tax_to_exclude} tax_filtered_table.qza ${asv_seq} tax_filtered_seq.qza tax_filtered_table.qzv ${metadata} tax_filtered_output tax_filtered_seq.qzv ${asv_tax_tsv} tax_filtered_table_with_tax.biom tax_filtered_table_with_tax.tsv completecmd &> q2_filtering_tax.log 2>&1
             """
        }
    }

    microDecon_table = params.filtering_tax_enable ? tax_filtered_table_tsv : biom_tsv_decontam

    /*
    * STEP 9 -  Decontaminate samples with MicroDecon
    */
    if (params.microDecon_enable) {
        process microDecon_step1 {
    
        	label 'microdecon_env'
    
        	publishDir "${params.outdir}/${params.microDecon_dirname}", mode: 'copy', pattern: 'decontaminated_ASV_table.tsv'
        	publishDir "${params.outdir}/${params.microDecon_dirname}", mode: 'copy', pattern: 'abundance_removed.txt'
        	publishDir "${params.outdir}/${params.microDecon_dirname}", mode: 'copy', pattern: 'ASV_removed.txt'
        	publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern: 'decontaminated_ASV_table.tsv'
        	publishDir "${params.outdir}/${params.report_dirname}/microDecon", mode: 'copy', pattern: 'decontaminated_ASV_table.tsv'
        	publishDir "${params.outdir}/${params.report_dirname}/microDecon", mode: 'copy', pattern: 'abundance_removed.txt'
        	publishDir "${params.outdir}/${params.report_dirname}/microDecon", mode: 'copy', pattern: 'ASV_removed.txt'
        	publishDir "${params.outdir}/${params.report_dirname}/version", mode: 'copy', pattern: 'v_*.txt'
        	publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'completecmd', saveAs : { complete_cmd_microDecon -> "cmd/${task.process}_complete.sh" }
    
        	input :
        		file microDecon_table from microDecon_table
    
        	output :
        		file 'decontaminated_ASV_table.tsv' into decontam_table, decontam_table_step2, decontam_table_step3
        		file 'abundance_removed.txt' into abund_removed
        		file 'ASV_removed.txt' into ASV_removed
        		file 'completecmd' into complete_cmd_microDecon
                        file 'v_microdecon.txt' into microdecon_version
    
        	when :
        		!params.stats_only && !params.dada2merge && params.microDecon_enable
    
        	shell :
        	"""
        	sed '1d' ${microDecon_table} > microDecon_table
        	sed -i 's/#OTU ID/ASV_ID/g' microDecon_table
        	microDecon.R microDecon_table ${params.control_list} ${params.nb_controls} ${params.nb_samples} decontaminated_ASV_table.tsv abundance_removed.txt ASV_removed.txt &> microDecon.log 2>&1
        	cp ${baseDir}/bin/microDecon.R completecmd &>> microDecon.log 2>&1
            Rscript -e "write(x=as.character(packageVersion('microDecon')), file='v_microdecon.txt')"
    
        	"""
        }
    
        process microDecon_step2 {
    
        	label 'qiime2_env'
    
        	publishDir "${params.outdir}/${params.microDecon_dirname}", mode: 'copy', pattern: 'decontaminated_ASV_table.qza'
        	publishDir "${params.outdir}/${params.report_dirname}/microDecon", mode: 'copy', pattern: 'decontaminated_ASV_table.qza'
    
        	input :
        		file table4microDecon from decontam_table_step2
    
        	output :
        		file 'decontaminated_ASV_table.qza' into decontam_table_qza, decontam_table_picrust2, decontam_table_ancom
    
        	when :
        		!params.stats_only && !params.dada2merge && params.microDecon_enable
    
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
        		file dada2_output from decontam_output
    
        	output :
        		file 'decontaminated_ASV_ID.txt' into decontam_ASV_ID
        		file 'decontaminated_ASV.fasta' into decontam_ASV_fasta
    
        	when :
        		!params.stats_only && !params.dada2merge && params.microDecon_enable
    
        	shell :
        	"""
        	cut -d \$'\t' -f1 ${decontam_table} | sed '1d' > decontaminated_ASV_ID.txt
        	seqtk subseq ${dada2_output}/sequences.fasta decontaminated_ASV_ID.txt > decontaminated_ASV.fasta
        	"""
        }
    
        /* Reformat decontaminated sequence to QIIME2 format */
        process microDecon_step4 {
    
        	label 'qiime2_env'
    
        	publishDir "${params.outdir}/${params.phylogeny_dirname}", mode: 'copy', pattern: '*.qza'
    
        	input :
        		file ASV_fasta from decontam_ASV_fasta
    
        	output :
        		file 'decontam_seqs.qza' into decontam_seqs_qza, decontam_seqs_phylo, decontam_seqs_picrust2, decontam_seqs_ancom
    
        	when :
        		!params.stats_only && !params.dada2merge && params.microDecon_enable
    
        	shell :
        	"""
                qiime tools import --input-path ${ASV_fasta} --output-path decontam_seqs.qza --type 'FeatureData[Sequence]'
        	"""
        }
    }
} else {
   
   if (workflow.profile.contains('test')) {
     longreadsmanifest = testmanifest.splitCsv(header: true, sep:'\t')
                                     .map { row -> tuple( row."sample-id", file(row."absolute-filepath")) }
     longreadstofasta = testmanifest.splitCsv(header: true, sep:'\t')
                                    .map { row -> file(row."absolute-filepath") }
   } else {
     longreadsmanifest = manifest.splitCsv(header: true, sep:'\t')
                                 .map { row -> tuple( row."sample-id", file(row."absolute-filepath")) }
     longreadstofasta = manifest.splitCsv(header: true, sep:'\t')
                                .map { row -> file(row."absolute-filepath") }
   }

   /* _______________________________________________________________ */
   /*                      Long reads part                            */
   /* _______________________________________________________________ */
   
   process lr_mapping {
     label 'lr_mapping_env'
   
     publishDir "${params.outdir}/${params.lr_mapping_dirname}", mode: 'copy', pattern: '*.bam'
     publishDir "${params.outdir}/${params.report_dirname}/version", mode: 'copy', pattern: 'v_*.txt'
   
     input:
       set sample, file(fastq) from longreadsmanifest
       file ready from ready_lr
   
     output:
       file "*.bam" into lr_mapped
       file 'lr_mapping.ok' into process_lr_mapping_report
       file 'lr_samtools-sort.log' into process_lr_sort_report
       file 'v_minimap2.txt' into v_minimap2_version
   
     shell:
       """
       minimap2 -t ${task.cpus} -K 25M -ax ${params.lr_type} -L ${params.lr_tax_fna} ${fastq} | samtools view -h -F0xe00 | samtools sort -o ${sample}.bam -O bam - &> lr_minimap2.log 2>&1
       touch lr_mapping.ok
       samtools index ${sample}.bam &> lr_samtools-sort.log 2>&1
       touch lr_samtools-sort.ok
       minimap2 --version > v_minimap2.txt
       """
   }

   process lr_getfasta {
      label 'seqtk_env'
      input : 
         file(fastq) from longreadstofasta

      output :
         file 'lr_sequences.fasta' into lr_sequences

      shell :
      """
         seqtk seq -a ${fastq} > 'lr_sequences.fasta'
      """
   }
 
   lr_sequences.collectFile(name : 'lr_sequences.fasta', newLine : false, storeDir : "${params.outdir}/${params.lr_mapping_dirname}")
               .subscribe {       println "Fasta sequences are saved to file : $it"       }

   process lr_get_taxonomy {
     label 'biopython_env'
 
     publishDir "${params.outdir}/${params.lr_taxonomy_dirname}", mode: 'copy', pattern: '*.tsv'
   
     input:
       file '*' from lr_mapped.collect()
   
     output:
       file 'samples.tsv' into lr_biom_tsv
       file 'lr_count_table.ok' into process_lr_taxonomy_report
   
     shell:
       """
       lr_count_table_minimap2.py -p "." -t "${params.lr_taxo_flat}" -r ${params.lr_rank} -o samples.tsv &> lr_count_table.log 2>&1
       touch lr_count_table.ok
       """
   }
}

if (!params.longreads) {
   seqs_phyloA = params.dada2merge ? merge_seqs_phylo : dada2_seqs_phylo
   seqs_phyloB = params.dbotu3_enable ? dbotu3_seqs_phylo : seqs_phyloA
   seqs_phyloC = params.filtering_tax_enable ? tax_filtered_seq : seqs_phyloB
   seqs_phylo = params.microDecon_enable ? decontam_seqs_phylo : seqs_phyloC

    /*
     * STEP 10 -  Run phylogeny construction
     */
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
    	    !params.stats_only
    
    	script :
    	"""
    	q2_phylogeny.sh ${repseqs_phylo} aligned_repseq.qza masked-aligned_repseq.qza tree.qza tree.log rooted_tree.qza tree_export_dir tree_export.log completecmd ${task.cpus} &>> q2_phylogeny.log 2>&1
    	cp tree_export_dir/tree.nwk tree.nwk &>> q2_phylogeny.log 2>&1
    	"""
    }

   table_picrust2A = params.dada2merge ? merge_table_picrust2 : dada2_table_picrust2
   table_picrust2B = params.dbotu3_enable ? dbotu3_table_picrust2 : table_picrust2A
   table_picrust2C = params.filtering_tax_enable ? tax_filtered_table_picrust2 : table_picrust2B
   table_picrust2 = params.microDecon_enable ? decontam_table_picrust2 : table_picrust2C
   
   seqs_picrust2A = params.dada2merge ? merge_seqs_picrust2 : dada2_seqs_picrust2
   seqs_picrust2B = params.dbotu3_enable ? dbotu3_seqs_picrust2 : seqs_picrust2A
   seqs_picrust2C = params.filtering_tax_enable ? tax_filtered_seq_picrust2 : seqs_picrust2B
   seqs_picrust2 = params.microDecon_enable ? decontam_seqs_picrust2 : seqs_picrust2C
   
   Channel
     .from(params.picrust_var)
     .splitCsv(sep : ',', strip : true)
     .flatten()
     .set { var_picrust2 }
   
   /*
    * STEP 11 -  Run functional predictions
    */
   if (params.picrust2_enable) {
   
       process q2_picrust2_analysis {
   
       	label 'qiime2_env'
   
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
       		!params.stats_only
   
       	script :
       	"""
       	q2_picrust2.sh ${table_picrust2} ${seqs_picrust2} q2-picrust2_output ${task.cpus} ${params.method} ${params.nsti} complete_picrust2_cmd &> q2_picrust2.log 2>&1
       	"""
       }
   
       /* Statistical analysis of functional predictions  */
       process q2_picrust2_stats {
   
           tag "$picrust_vars"
           label 'r_stats_env'
           label 'internet_access'
   
           publishDir "${params.outdir}/${params.report_dirname}/picrust2_output", mode: 'copy', pattern: '*functional_predictions_NMDS*'
           publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'complete_picrust2_stats_cmd', saveAs : { complete_picrust2_stats_cmd -> "cmd/${task.process}_complete.sh" }
   
           input :
               file ec_metagenome from EC_predictions_tsv
               file ko_metagenome from KO_predictions_tsv
               file metacyc_predictions_ from pathway_predictions_tsv
               file metadata from metadata4picrust2
               each picrust_vars from var_picrust2
   
           output :
               file '*functional_predictions_NMDS*' into functional_pred_NMDS
               file 'complete_picrust2_stats_cmd' into complete_picrust2_stats_cmd
   
           when :
               !params.stats_only
   
           script :
           """
           functional_predictions.R ec_metagenome_predictions_with-descriptions.tsv ko_metagenome_predictions_with-descriptions.tsv pathway_abundance_predictions_with-descriptions.tsv ${metadata} ${picrust_vars} functional_predictions_NMDS_${picrust_vars} ${params.microDecon_enable} ${params.control_list} &> picrust2_stats.log 2>&1
           cp ${baseDir}/bin/functional_predictions.R complete_picrust2_stats_cmd &>> picrust2_stats.log 2>&1
           """
       }
   }
   
   table_ancomA = params.dada2merge ? merge_table_ancom : dada2_table_ancom
   table_ancomB = params.dbotu3_enable ? dbotu3_table_ancom : table_ancomA
   table_ancomC = params.filtering_tax_enable ? tax_filtered_table_ancom : table_ancomB
   table_ancom = params.microDecon_enable ? decontam_table_ancom : table_ancomC
   
   Channel
     .from(params.ancom_var)
     .splitCsv(sep : ',', strip : true)
     .flatten()
     .set { ancom_var_list }
   
   /*
    * STEP 12 -  Differential abundance testing with ANCOM
    */
   process q2_ancom {
   
       tag "$ancom_var"
       label 'qiime2_env'
   
       publishDir "${params.outdir}/${params.ancom_dirname}", mode: 'copy', pattern: '*.qz*'
       publishDir "${params.outdir}/${params.report_dirname}/ancom_output", mode: 'copy', pattern: 'export_ancom_*'
       publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'completecmd_ancom', saveAs : { completecmd_ancom -> "cmd/${task.process}_complete.sh" }
   
       input :
           file table4ancom from table_ancom
           file metadata from metadata4ancom
           file taxonomy4ancom from data_taxonomy_ancom
           each ancom_var from ancom_var_list
   
       output :
           file 'compo_table*.qza' into compo_table
           file 'ancom_*.qzv' into ancom_table
           file 'export_ancom_*' into export_ancom
           file 'collapsed_table_*.qza' into collapsed_taxolevel_table
           file 'completecmd_ancom' into completecmd_ancom, completcmd_ancom4compress
   
       when :
           !params.stats_only && params.ancom_enable
   
       script :
       """
       q2_ANCOM.sh ${table4ancom} compo_table.qza ${metadata} ${ancom_var} ancom_${ancom_var}.qzv export_ancom_${ancom_var} ${taxonomy4ancom} collapsed_table_family.qza compo_table_family.qza ancom_${ancom_var}_family.qzv export_ancom_${ancom_var}_family collapsed_table_genus.qza compo_table_genus.qza ancom_${ancom_var}_genus.qzv export_ancom_${ancom_var}_genus completecmd_ancom &> q2_ancom.log 2>&1
       """
   }
}

if (params.longreads) {
   tsv = lr_biom_tsv
   newick = "none"
} else {   
   tsvA = params.filtering_tax_enable ? tax_filtered_table_tsv_stats : biom_tsv
   tsvB = params.microDecon_enable ? decontam_table : tsvA
   tsv = params.stats_only ? tsv_only : tsvB
   newick = params.stats_only ? newick_only : newick_phylo
}

/*
 * STEP 13 -  Prepare data for statistics steps
 */
process prepare_data_for_stats {

    label 'r_stats_env'

    publishDir "${params.outdir}/${params.report_dirname}/R/DATA", mode: 'copy', pattern : '*.tsv'
    publishDir "${params.outdir}/${params.report_dirname}/R/DATA", mode: 'copy', pattern : '*.rds'
    publishDir "${params.outdir}/${params.report_dirname}/version", mode: 'copy', pattern : 'v_*.txt'

    input :
        file metadata from metadata4stats
        file biom_tsv from tsv
        file newick_tree from newick

    output :
        file 'table_with_taxo_for_stats.tsv' into biom_tsv_stats
        file 'metadata_stats.tsv' into metadata_stats, metadata_beta, metadata_beta_rarefied, metadata_beta_deseq2, metadata_beta_css
        file 'phyloseq.rds' into phyloseq_rds, phyloseq_rds_alpha, phyloseq_rds_beta, phyloseq_rds_beta_rarefied, phyloseq_rds_beta_deseq2, phyloseq_rds_beta_css,phyloseq_rds_set
        file 'version_ok' into version_collected
        file 'v_*.txt' into r_lib_version
 
    script :
    """
    if [ ${params.longreads} = "true" ];
    then
        Rscript --vanilla ${baseDir}/bin/create_phyloseq_obj_longreads.R phyloseq.rds ${biom_tsv} ${metadata} ${params.lr_rank} ${params.kingdom} table_with_taxo_for_stats.tsv &>> stats_prepare_data_lr.log 2&>1
    else
        prepare_data_for_stats.sh ${metadata} ${biom_tsv} table_with_taxo_for_stats.tsv metadata_stats.tsv ${params.microDecon_enable} ${params.stats_only} &> stats_prepare_data.log 2&>1
        Rscript --vanilla ${baseDir}/bin/create_phyloseq_obj.R phyloseq.rds table_with_taxo_for_stats.tsv metadata_stats.tsv ${params.microDecon_enable} ${params.control_list} ${newick_tree} ${params.remove_sample} ${params.sample_to_remove} &>> stats_prepare_data.log 2&>1
    fi
    ## get statistics libraries version for report
    Rscript -e "write(x=as.character(paste0(R.Version()[c('major','minor')], collapse = '.')), file='v_R.txt')"
    Rscript -e "library(dplyr); write(x=as.character(packageVersion('dplyr')), file='v_dplyr.txt')"
    Rscript -e "library(stringr); write(x=as.character(packageVersion('stringr')), file='v_stringr.txt')"
    Rscript -e "library(phyloseq); x=as.character(packageVersion('phyloseq')); write(x, file='v_phyloseq.txt')"
    Rscript -e "library(phangorn); write(x=as.character(packageVersion('phangorn')), file='v_phangorn.txt')"
    Rscript -e "library(ggplot2); write(x=as.character(packageVersion('ggplot2')), file='v_ggplot2.txt')"
    Rscript -e "library(svglite); write(x=as.character(packageVersion('svglite')), file='v_svglite.txt')"
    Rscript -e "library(vegan); write(x=as.character(packageVersion('vegan')), file='v_vegan.txt')"
    Rscript -e "library(RColorBrewer); write(x=as.character(packageVersion('RColorBrewer')), file='v_RColorBrewer.txt')"
    Rscript -e "library(tidyr); write(x=as.character(packageVersion('tidyr')), file='v_tidyr.txt')"
    Rscript -e "library(gridExtra); write(x=as.character(packageVersion('gridExtra')), file='v_gridExtra.txt')"
    Rscript -e "library(egg); write(x=as.character(packageVersion('egg')), file='v_egg.txt')"
    Rscript -e "library(reshape2); write(x=as.character(packageVersion('reshape2')), file='v_reshape2.txt')"
    Rscript -e "library(BiocManager); write(x=as.character(packageVersion('BiocManager')), file='v_BiocManager.txt')"
    Rscript -e "library(microbiome); write(x=as.character(packageVersion('microbiome')), file='v_microbiome.txt')"
    Rscript -e "library(dendextend); write(x=as.character(packageVersion('dendextend')), file='v_dendextend.txt')"
    Rscript -e "library(metagenomeSeq); write(x=as.character(packageVersion('metagenomeSeq')), file='v_metagenomeSeq.txt')"
    Rscript -e "library(DESeq2); write(x=as.character(packageVersion('DESeq2')), file='v_DESeq2.txt')"
    Rscript -e "library(UpSetR); write(x=as.character(packageVersion('UpSetR')), file='v_UpSetR.txt')"
    touch 'version_ok'
    """
}

Channel
  .from(params.alpha_div_group)
  .splitCsv(sep : ',', strip : true)
  .flatten()
  .set { alpha_list_var }

/*
 * STEP 14 -  Alpha diversity community statistics analysis
 */
if (params.stats_alpha_enable) {
    process stats_alpha {
    
        tag "$alpha_var"
        label 'r_stats_env'
        label 'internet_access'
    
        publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/alpha_diversity", mode: 'copy', pattern : 'alpha_div_values.txt'
        publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/alpha_diversity", mode: 'copy', pattern : 'index_significance_tests.txt'
        publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/alpha_diversity/diversity_index", mode: 'copy', pattern : 'alpha_div_plots*'
        publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/alpha_diversity/diversity_barplots/${alpha_var}", mode: 'copy', pattern : 'barplot_*'
        publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/alpha_diversity", mode: 'copy', pattern : 'rarefaction_curve*'
    
        input :
           file phyloseq_rds from phyloseq_rds_alpha
           each alpha_var from alpha_list_var
    
        output :
            file 'alpha_div_values.txt' into alpha_div_values
            file "alpha_div_plots_${alpha_var}*" into alpha_div_plots
            file 'index_significance_tests.txt' into index_significance_tests
            file "barplot_*_${alpha_var}*" into barplots
            file 'rarefaction_curve*' into rarefaction_curve
            file 'process_alpha_report.ok' into process_alpha_report
    
        when :
            params.stats_alpha_enable
    
        shell :
        """
        Rscript --vanilla ${baseDir}/bin/alpha_diversity.R ${phyloseq_rds} alpha_div_values.txt alpha_div_plots_${alpha_var} ${params.kingdom} ${params.taxa_nb} barplot_phylum_${alpha_var} barplot_class_${alpha_var} barplot_order_${alpha_var} barplot_family_${alpha_var} barplot_genus_${alpha_var} ${alpha_var} index_significance_tests.txt $workflow.projectDir rarefaction_curve ${params.longreads} ${params.lr_rank} &> stats_alpha_diversity.log 2>&1
        touch process_alpha_report.ok
        """
    }
}

Channel
  .from(params.beta_div_var)
  .splitCsv(sep : ',', strip : true)
  .flatten()
  .into { beta_var_nn ; beta_var_rare ; beta_var_deseq ; beta_var_css }

/*
 * STEP 15 -  Beta diversity (non normalized) community statistics analysis
 */
if (params.stats_beta_enable) {
    process stats_beta {
    
        tag "$beta_var"
        label 'r_stats_env'
        label 'internet_access'
    
        publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_non_normalized/NMDS", mode: 'copy', pattern : 'NMDS*'
        publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_non_normalized/PCoA", mode: 'copy', pattern : 'PCoA*'
        publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_non_normalized/hclustering", mode: 'copy', pattern : 'hclustering*'
        publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_non_normalized", mode:'copy', pattern : 'variance_significance_tests_*'
        publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_non_normalized/ExpVar", mode:'copy', pattern : 'pie_ExpVar_*'
        publishDir "${params.outdir}/${params.report_dirname}/R/DATA", mode: 'copy', pattern : '*.tsv'
    
        input :
            file phyloseq_rds from phyloseq_rds_beta
            file metadata from metadata_beta
            each beta_var from beta_var_nn
    
        output :
            file "NMDS_${beta_var}*" into NMDS
            file "PCoA_${beta_var}*" into PCoA
            file "hclustering_${beta_var}*" into hclustering
            file 'variance_significance_tests_*' into variance_significance_tests
            file 'pie_ExpVar_*' into pie_ExpVar
            file 'process_beta_report.ok' into process_beta_report
    
        when :
            params.stats_beta_enable
    
        shell :
        """
        Rscript --vanilla ${baseDir}/bin/beta_diversity.R ${phyloseq_rds} ${beta_var} ${metadata} $workflow.projectDir NMDS_${beta_var}_ PCoA_${beta_var}_ ${params.hc_method} hclustering_${beta_var}_ variance_significance_tests_ pie_ExpVar_ ${params.longreads} &> stats_beta_diversity.log 2>&1
        touch process_beta_report.ok
        """
    }
    
    /*
     * STEP 16 -  Beta diversity (rarefied) community statistics analysis
     */
    process stats_beta_rarefied {
    
        tag "$beta_var"
        label 'r_stats_env'
        label 'internet_access'
    
        publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_rarefied/NMDS", mode: 'copy', pattern : 'NMDS*'
        publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_rarefied/PCoA", mode: 'copy', pattern : 'PCoA*'
        publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_rarefied/hclustering", mode: 'copy', pattern : 'hclustering*'
        publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_rarefied", mode: 'copy', pattern : 'variance_significance_tests_rarefied_*'
        publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_rarefied/ExpVar", mode:'copy', pattern : 'pie_ExpVar_rarefied_*'
        publishDir "${params.outdir}/${params.report_dirname}/R/DATA", mode: 'copy', pattern : '*.tsv'
    
        input :
            file phyloseq_rds from phyloseq_rds_beta_rarefied
            file metadata from metadata_beta_rarefied
            each beta_var from beta_var_rare
    
        output :
            file 'Final_rarefied_ASV_table_with_taxonomy.tsv' into final_rarefied_ASV_table_with_taxonomy
            file "NMDS_rarefied_${beta_var}*" into NMDS_rarefied
            file "PCoA_rarefied_${beta_var}*" into PCoA_rarefied
            file "hclustering_rarefied_${beta_var}*" into hclustering_rarefied
            file 'variance_significance_tests_rarefied_*' into variance_significance_tests_rarefied
            file 'pie_ExpVar_rarefied_*' into pie_ExpVar_rarefied
            file 'process_beta_report_rarefied.ok' into process_beta_report_rarefied
    
        when :
            params.stats_beta_enable
    
        shell :
        """
        Rscript --vanilla ${baseDir}/bin/beta_diversity_rarefied.R ${phyloseq_rds} Final_rarefied_ASV_table_with_taxonomy.tsv ${beta_var} ${metadata} $workflow.projectDir NMDS_rarefied_${beta_var}_ PCoA_rarefied_${beta_var}_ ${params.hc_method} hclustering_rarefied_${beta_var}_ variance_significance_tests_rarefied_ pie_ExpVar_rarefied_ ${params.longreads} &> stats_beta_diversity_rarefied.log 2>&1
        touch process_beta_report_rarefied.ok
        """
    }
    
    /*
     * STEP 17 -  Beta diversity (DESeq2) community statistics analysis
     */
    process stats_beta_deseq2 {
    
        tag "$beta_var"
        label 'r_stats_env'
        label 'internet_access'
    
        publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_DESeq2/NMDS", mode: 'copy', pattern : 'NMDS*'
        publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_DESeq2/PCoA", mode: 'copy', pattern : 'PCoA*'
        publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_DESeq2/hclustering", mode: 'copy', pattern : 'hclustering*'
        publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_DESeq2", mode: 'copy',pattern : 'variance_significance_tests_DESeq2_*'
        publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_DESeq2/ExpVar", mode: 'copy',pattern : 'pie_ExpVar_DESeq2_*'
        publishDir "${params.outdir}/${params.report_dirname}/R/DATA", mode: 'copy', pattern : '*.tsv'
    
        input :
            file phyloseq_rds from phyloseq_rds_beta_deseq2
            file metadata from metadata_beta_deseq2
            each beta_var from beta_var_deseq
    
        output :
            file 'Final_DESeq2_ASV_table_with_taxonomy.tsv' into final_deseq2_ASV_table_with_taxonomy
            file "NMDS_DESeq2_${beta_var}*" into NMDS_deseq2
            file "PCoA_DESeq2_${beta_var}*" into PCoA_deseq2
            file "hclustering_DESeq2_${beta_var}*" into hclustering_deseq2
            file 'variance_significance_tests_DESeq2_*' into variance_significance_tests_DESeq2
            file 'pie_ExpVar_DESeq2_*' into pie_ExpVar_DESeq2
            file 'process_beta_report_DESeq2.ok' into process_beta_report_DESeq2
    
        when :
            params.stats_beta_enable
    
        shell :
        """
        Rscript --vanilla ${baseDir}/bin/beta_diversity_deseq2.R ${phyloseq_rds} Final_DESeq2_ASV_table_with_taxonomy.tsv ${beta_var} ${metadata} $workflow.projectDir NMDS_DESeq2_${beta_var}_ PCoA_DESeq2_${beta_var}_ ${params.hc_method} hclustering_DESeq2_${beta_var}_ variance_significance_tests_DESeq2_ pie_ExpVar_DESeq2_ ${params.longreads} &> stats_beta_diversity_deseq2.log 2>&1
        touch process_beta_report_DESeq2.ok
        """
    }
    
    /*
     * STEP 18 -  Beta diversity (CSS) community statistics analysis
     */
    process stats_beta_css {
    
        tag "$beta_var"
        label 'r_stats_env'
        label 'internet_access'
    
        publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_CSS/NMDS", mode: 'copy', pattern : 'NMDS*'
        publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_CSS/PCoA", mode: 'copy', pattern : 'PCoA*'
        publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_CSS/hclustering", mode: 'copy', pattern : 'hclustering*'
        publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_CSS", mode: 'copy', pattern : 'variance_significance_tests_CSS_*'
        publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/beta_diversity_CSS/ExpVar", mode: 'copy', pattern : 'pie_ExpVar_CSS_*'
        publishDir "${params.outdir}/${params.report_dirname}/R/DATA", mode: 'copy', pattern : '*.tsv'
    
        input :
            file phyloseq_rds from phyloseq_rds_beta_css
            file metadata from metadata_beta_css
            each beta_var from beta_var_css
    
        output :
            file 'Final_CSS_ASV_table_with_taxonomy.tsv' into final_css_ASV_table_with_taxonomy
            file "NMDS_CSS_${beta_var}*" into NMDS_css
            file "PCoA_CSS_${beta_var}*" into PCoA_css
            file "hclustering_CSS_${beta_var}*" into hclustering_css
            file 'variance_significance_tests_CSS_*' into variance_significance_tests_CSS
            file 'pie_ExpVar_CSS_*' into pie_ExpVar_CSS
            file 'process_beta_report_CSS.ok' into process_beta_report_CSS
    
        when :
            params.stats_beta_enable
    
        shell :
        """
        Rscript --vanilla ${baseDir}/bin/beta_diversity_css.R ${phyloseq_rds} Final_CSS_ASV_table_with_taxonomy.tsv ${beta_var} ${metadata} $workflow.projectDir NMDS_CSS_${beta_var}_ PCoA_CSS_${beta_var}_ ${params.hc_method} hclustering_CSS_${beta_var}_ variance_significance_tests_CSS_ pie_ExpVar_CSS_ ${params.longreads} &> stats_beta_diversity_css.log 2>&1
        touch process_beta_report_CSS.ok
        """
    }
}

Channel
  .from(params.desc_comp_crit)
  .splitCsv(sep : ',', strip : true)
  .flatten()
  .set { desc_comp_list }

/*
 * STEP 19 -  Descriptive comparisons statistics analysis
 */
if (params.stats_desc_comp_enable) {
    process stats_desc_comp {
    
        tag "$desc_comp_var"
        label 'r_stats_env'
    
        publishDir "${params.outdir}/${params.report_dirname}/R/FIGURES/descriptive_comparison", mode: 'copy', pattern : 'upset_plot*'
    
        input :
            file phyloseq_rds from phyloseq_rds_set
            each desc_comp_var from desc_comp_list
    
        output :
            file "upset_plot_${desc_comp_var}*" into upset_plot
            file 'process_desc_comp_report.ok' into process_desc_comp_report
    
        when :
            params.stats_desc_comp_enable
    
        shell :
        """
        Rscript --vanilla ${baseDir}/bin/desc_comp.R ${phyloseq_rds} ${desc_comp_var} upset_plot_${desc_comp_var} ${params.desc_comp_tax_level} &> stats_desc_comp.log 2>&1
        touch process_desc_comp_report.ok
        """
    }
}

SAMBAtemplate_ch = params.report_enable ? Channel.fromPath(params.SAMBAtemplate, checkIfExists:true) : Channel.empty()
SAMBAcss_ch = params.report_enable ? Channel.fromPath(params.SAMBAcss, checkIfExists:true) : Channel.empty()
SAMBAlogo_ch = params.report_enable ? Channel.fromPath(params.SAMBAlogo, checkIfExists:true) : Channel.empty()
SAMBAwf_ch = params.report_enable ? Channel.fromPath(params.SAMBAwf, checkIfExists:true) : Channel.empty()
betastats_reportok = params.stats_beta_enable ? process_beta_report.concat( process_beta_report_CSS, process_beta_report_DESeq2, process_beta_report_rarefied ) : Channel.empty()
SAMBAreport_okstats_alpha = params.stats_alpha_enable ? process_alpha_report : Channel.from('report_without_stats_alpha_ok')
SAMBAreport_okstats_beta = params.stats_beta_enable ? betastats_reportok : Channel.from('report_without_stats_beta_ok')
SAMBAreport_okdesc_comp = params.stats_desc_comp_enable ? process_desc_comp_report : Channel.from('report_without_desc_comp_ok')

SAMBAreport_okpicrust2 = params.picrust2_enable ? complete_picrust2_stats_cmd : Channel.from('report_without_picrust2_ok')

SAMBAreport_okancom = params.ancom_enable ? completecmd_ancom : Channel.from('report_without_ancom_ok')

/*
 * STEP 20 -  Save user parameters of the workflow and Generate analysis report
 */
if (params.report_enable) {
    process report {

        label 'jinja2_env'

        publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'SAMBA_report.html'
        publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'style.css'
        publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'nfcore-samba_logo.png'
        publishDir "${params.outdir}/${params.report_dirname}", mode: 'copy', pattern : 'samba_wf.png'
        publishDir "${params.outdir}/conf", mode: 'copy', pattern : 'data.json'

        input :
            file SAMBAtemplate from SAMBAtemplate_ch
            file SAMBAcss from SAMBAcss_ch
            file SAMBAreport_okstats_alpha from SAMBAreport_okstats_alpha
            file SAMBAreport_okstats_beta from SAMBAreport_okstats_beta
            file SAMBAreport_okdesc_comp from SAMBAreport_okdesc_comp
            file SAMBAreport_okpicrust2 from SAMBAreport_okpicrust2
            file SAMBAreport_okancom from SAMBAreport_okancom
            file logo from SAMBAlogo_ch
            file wf_image from SAMBAwf_ch
            file 'version_ok' from version_collected

       output :
            file 'style.css' into SAMBA_css_output
            file 'SAMBA_report.html' into Report
            file 'nfcore-samba_logo.png' into SAMBAlogo_output
            file 'samba_wf.png' into wf_image_output
            file 'data.json' into data_json

        when :
           params.report_enable

        script :
        """
        #!/usr/bin/env python
        import os, json
        from shutil import copyfile

        data = {}
        data["projectName"] = '$params.projectName'
        data["run_date"] = '$date'
        data["singleEnd"] = '$params.singleEnd'
        data["manifest"] = '$params.input_manifest'
        data["metadata"] = '$params.input_metadata'
        data["outdir"] = '$params.outdir'
        data["steps"] = {}
        data["steps"]["data_integrity_enable"] = '$params.data_integrity_enable'
        data["steps"]["dbotu3_enable"] = '$params.dbotu3_enable'
        data["steps"]["filtering_tax_enable"] = '$params.filtering_tax_enable'
        data["steps"]["microDecon_enable"] = '$params.microDecon_enable'
        data["steps"]["picrust2_enable"] = '$params.picrust2_enable'
        data["steps"]["stats_alpha_enable"] = '$params.stats_alpha_enable'
        data["steps"]["stats_beta_enable"] = '$params.stats_beta_enable'
        data["steps"]["stats_desc_comp_enable"] = '$params.stats_desc_comp_enable'
        data["steps"]["report_enable"] = '$params.report_enable'
        data["steps"]["stats_only"] = '$params.stats_only'
        data["steps"]["dada2merge"] = '$params.dada2merge'
        data["steps"]["longreads"] = '$params.longreads'

        data["integrity"] = {}
        data["integrity"]["primer_filter"] = '$params.primer_filter' or None

        data["cutadapt"] = {}
        data["cutadapt"]["primerF"] = '$params.primerF' or None
        data["cutadapt"]["primerR"] = '$params.primerR' or None
        data["cutadapt"]["errorRate"] = '$params.errorRate' or None
        data["cutadapt"]["overlap"] = '$params.overlap' or None

        data["dada2"] = {}
        data["dada2"]["FtrimLeft"] = '$params.FtrimLeft' or None
        data["dada2"]["RtrimLeft"] = '$params.RtrimLeft' or None
        data["dada2"]["FtruncLen"] = '$params.FtruncLen' or None
        data["dada2"]["RtruncLen"] = '$params.RtruncLen' or None
        data["dada2"]["FmaxEE"] = '$params.FmaxEE' or None
        data["dada2"]["RmaxEE"] = '$params.RmaxEE' or None
        data["dada2"]["minQ"] = '$params.minQ' or None
        data["dada2"]["chimeras"] = '$params.chimeras' or None

        data["dada2merge"] = {}
        data["dada2merge"]["merge_tabledir"] = '$params.merge_tabledir' or None
        data["dada2merge"]["merge_repseqsdir"] = '$params.merge_repseqsdir' or None

        data["dbotu3"] = {}
        data["dbotu3"]["gen_crit"] = '$params.gen_crit' or None
        data["dbotu3"]["abund_crit"] = '$params.abund_crit' or None
        data["dbotu3"]["pval_crit"] = '$params.pval_crit' or None

        data["taxonomy"] = {}
        data["taxonomy"]["database"] = '$params.database' or None
        data["taxonomy"]["seqs_db"] = '$params.seqs_db' or None
        data["taxonomy"]["taxo_db"] = '$params.taxo_db' or None
        data["taxonomy"]["extract_db"] = '$params.extract_db'
        data["taxonomy"]["confidence"] = '$params.confidence' or None

        data["picrust2"] = {}
        data["picrust2"]["method"] = '$params.method' or None
        data["picrust2"]["nsti"] = '$params.nsti' or None

        data["tax_filtering"] = {}
        data["tax_filtering"]["tax_to_exclude"] = '$params.tax_to_exclude' or None

        data["microdecon"] = {}
        data["microdecon"]["control_list"] = '$params.control_list' or None
        data["microdecon"]["nb_controls"] = '$params.nb_controls' or None
        data["microdecon"]["nb_samples"] = '$params.nb_samples' or None

        data["lr"] = {}
        data["lr"]["lr_type"] = '$params.lr_type' or None
        data["lr"]["lr_rank"] = '$params.lr_rank' or None
        data["lr"]["lr_tax_fna"] = '$params.lr_tax_fna' or None
        data["lr"]["lr_tax_flat"] = '$params.lr_taxo_flat' or None

        data["stats"] = {}
        data["stats"]["ancom_var"] = '$params.ancom_var' or None
        data["stats"]["kingdom"] = '$params.kingdom' or None
        data["stats"]["taxa_nb"] = '$params.taxa_nb' or None
        data["stats"]["hc_method"] = '$params.hc_method' or None
        data["stats"]["alpha_div_group"] = '$params.alpha_div_group' or None
        data["stats"]["beta_div_var"] = '$params.beta_div_var' or None
        data["stats"]["desc_comp_crit"] = '$params.desc_comp_crit' or None
        data["stats"]["inasv_table"] = '$params.inasv_table' or None
        data["stats"]["innewick"] = '$params.innewick' or None

        #software versions
        data["soft"] = {}
        data["soft"]["samba"] = '$workflow.manifest.version'
        data["soft"]["nextflow"] = '$workflow.nextflow.version'
        for elt in os.listdir("${params.outdir}/${params.report_dirname}/version"):
            if elt.endswith(".txt"):
                f=open(os.path.join("${params.outdir}/${params.report_dirname}/version",elt),"r")
                tool_name=elt.replace("v_","").replace(".txt","")
                version=f.readline().rstrip()
                data["soft"][tool_name] = version

        with open('data.json', 'w') as outfile:
           json.dump(data, outfile, sort_keys=True, indent=4)

        os.system("${baseDir}/bin/SAMBAreport.py -t ${SAMBAtemplate} -p ${params.outdir}/${params.report_dirname} -c 'data.json'")
        os.system("cp ${wf_image} samba_wf.png")

        """
    }
}

/*
 * STEP 21 -  Compress final report directory
 */
if (params.compress_result) {
    process compress_result {
        
        input :
            file compress_ok from Report

        output :
            file "SAMBA_report.zip" into report_zip

        when :
           params.compress_result
 
        script :
        """
        cd ${params.outdir}/${params.report_dirname} && zip -r SAMBA_report.zip * 
        cd -
        ln -s ${params.outdir}/${params.report_dirname}/SAMBA_report.zip
        """
    }
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[samba workflow] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[samba workflow] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[samba workflow] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            [ 'mail', '-s', subject, email_address ].execute() << email_txt
            log.info "[samba workflow] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[samba workflow]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[samba workflow]${c_red} Pipeline completed with errors${c_reset}-"
    }
}

/* Other functions */
def SeBiMERHeader() {
    // Log colors ANSI codes
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    -${c_cyan}--------------------------------------------------${c_reset}-
    ${c_blue}    __  __  __  .       __  __  ${c_reset}
    ${c_blue}   \\   |_  |__) | |\\/| |_  |__)  ${c_reset}
    ${c_blue}  __\\  |__ |__) | |  | |__ |  \\  ${c_reset}
                                            ${c_reset}
    ${c_yellow}  samba v${workflow.manifest.version}: Standardized and Automated MetaBarcoding Analyses workflow${c_reset}
    -${c_cyan}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}
