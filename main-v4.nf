#!/usr/bin/env nextflow
/*
=======================================================================================
                                  SAMBA workflow                                  
=======================================================================================
 Standardized and Automated MetaBarcoding Analyses workflow
 #### Homepage / Documentation
 https://gitlab.ifremer.fr/bioinfo/SAMBA-nextflow
 https://github.com/ifremer-bioinformatics/samba
---------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

def helpMessage() {
    // Add to this help message with new command line parameters
    log.info SeBiMERHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline after filling the config file corresponding to your analysis as follows:

	nextflow run main.nf -profile <illumina_test/illumina/longreadstest>,singularity [-c <institute_config_file>]

	Mandatory:
	--excel_sample_file		[path]	Path to the XLS input file (EXCEL 97-2004) containing the manifest and metadata sheets.
	--longreads			[bool]	Set to true to specify that the inputs are long reads (Nanopore/Pacbio) (default = false for illumina short reads).
	--singleEnd			[bool]	Set to true to specify that the inputs are single-end reads (default = false).

	Other options:
	--outdir			[path]	The output directory where the results will be saved.
	-name 				[str]	Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
	--projectName			[str]	Name of the project.

	@@ PROCESS OPTIONS FOR ILLUMINA ANALYSIS @@

	Data integrity:
	--data_integrity_enable		[bool]	Data integrity checking step. Set to false to deactivate this step (default = true).
	--primer_filter			[str]	Percentage of primers supposed to be found in raw reads (default = 70).

	Cutadapt - primer removal:
	--cutadapt_enable		[bool]	Primer removal process. Set to false to deactivate this step (default = true).
	--primerF			[str]	Forward primer (to be used in Cutadapt cleaning step).
	--primerR			[str]	Reverse primer (to be used in Cutadapt cleaning step).
	--errorRate			[str]	Cutadapt error rate allowed to match primers (default = 0.1).

	FIGARO - optimizing trimming parameters for DADA2:
	--figaro_enable			[bool]	Process to identify optimal trimming parameters for DADA2. Set to false to deactivate this step (default = true).
	--raw_data_dir			[path]	Path to raw data directory.
	--amplicon_length		[str]	Length of expected amplicons.

	DADA2 - ASV inference:
	--FtrimLeft			[str]	The number of nucleotides to remove from the start of each forward read (default : 0 = no trimming).
	--RtrimLeft			[str]	The number of nucleotides to remove from the start of each reverse read (default : 0 = no trimming).
	--FtruncLen			[str]	Truncate forward reads after FtruncLen bases. Reads shorter than this are discarded (default : 0 = no trimming).
	--RtruncLen			[str]	Truncate reverse reads after RtruncLen bases. Reads shorter than this are discarded (default : 0 = no trimming).
	--FmaxEE			[str]	Forward reads with higher than max "expected errors" will be discarded (default = 2).
	--RmaxEE			[str]	Reverse with higher than max "expected errors" will be discarded (default = 2).
	--truncQ			[str]	Truncate reads at the first instance of a quality score less than or equal to minQ (default = 2).
	--pooling_method		[str]	Method used to pool samples for denoising. Default = "independant". Set to "pseudo" if you want to approximate pooling of samples (see DADA2 documentation).
	--chimeras_method		[str]	Chimera detection method : default = "consensus". Set to "pooled" if the samples in the sequence table are all pooled together for bimera identification (see DADA2 documentation).

	dbOTU3 - Distribution based-clustering:
	--dbotu3_enable			[bool]	Distribution based-clustering step. Set to false to deactivate
 this step (default = true).
	--gen_crit			[str]	dbOTU3 Genetic criterion (default = 0.1).
	--abund_crit			[str]	dbOTU3 Abundance criterion (default = 10).
	--pval_crit			[str]	dbOTU3 P-value criterion (default = 0.0005).

	Taxonomic assignation:
	--database			[file]	Path to a trained Naive Bayes QIIME 2 classifier.
	--confidence			[str]	Confidence threshold for limiting taxonomic depth. Set to "disable" to disable confidence calculation, or 0 to calculate confidence but not apply it to limit the taxonomic depth of the assignments (default = 0.7).

	Filter ASV table and ASV sequences based on taxonomy:
	--filter_table_by_tax_enable	[bool]	Set to true to filter ASV table and ASV sequences based on taxonomic assignation (default = false).
	--tax_to_exclude		[str]	List of taxa you want to exclude (comma-separated list).

	Filter ASV table and ASV sequences based on data:
	--filter_table_by_data_enable	[bool]	Set to true to filter ASV table and ASV sequences based on data (sample ID and/or frequency/contingency) (default = false).
	--filter_by_id			[bool]	Set to true if you want to filter data based on a list of sample IDs (default = false ; only if filter_table_by_data is enable).
	--list_sample_to_remove		[str]	List of sample IDs to remove (comma-separated list).
	--filter_by_frequency		[bool]	Set to true if you want to filter data based on frequency and contingency (default = false ; only if filter_table_by_data is enable).
	--min_frequency_sample		[int]	Minimum desired total sequence count within each sample (default = 1).
	--min_frequency_asv		[int]	Minimum desired abundance for each ASV (default = 1).
	--contingency_asv		[int]	Minimum number of samples in which an ASV must be found to be kept (default = 1).

	Filter contaminant ASVs:
	--filter_contaminants_enable	[bool]	Sample decontamination step. Set to true to activate this step (default = false)
	--list_control_samples		[str]	List of control sample IDs (comma-separated list).

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

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

//Copy config files to output directory for each run
base_params_file = file("${baseDir}/conf/base.config", checkIfExists: true)
base_params_file.copyTo("${params.outdir}/00_pipeline_config/base.config")
if (workflow.profile.contains('illumina')) {
    custom_params_file = file("${baseDir}/conf/illumina.config", checkIfExists: true)
    custom_params_file.copyTo("${params.outdir}/00_pipeline_config/illumina.config")
}
if (workflow.profile.contains('illumina_test')) {
    shortreadstest_params_file = file("${baseDir}/conf/illumina_test.config", checkIfExists: true)
    shortreadstest_params_file.copyTo("${params.outdir}/00_pipeline_config/illumina_test.config")
}
if (workflow.profile.contains('longreadstest')) {
   longreadstest_params_file = file("${baseDir}/conf/longreadstest.config", checkIfExists: true)
   longreadstest_params_file.copyTo("${params.outdir}/00_pipeline_config/longreadstest.config")
}

/*
 * PIPELINE INFO
 */

// Header log info
log.info SeBiMERHeader()
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name'] = custom_runName ?: workflow.runName
summary['Project Name'] = params.projectName
summary['User'] = workflow.userName
summary['Launch dir'] = workflow.launchDir
summary['Working dir'] = workflow.workDir
summary['Output dir'] = params.outdir
summary['Profile'] = workflow.profile
if (workflow.profile == 'illumina') summary['Data type'] = 'Illumina short reads'
if (workflow.profile == 'illumina_test') summary['Data type'] = 'Illumina short reads test workflow'
if (workflow.profile == 'longreadstest') summary['Data type'] = 'Long reads analysis'
summary['Sample Input Excel File'] = params.excel_sample_file
if (params.data_integrity_enable) summary['Data integrity'] = "Data integrity checking process enabled"
if (params.cutadapt_enable) summary['Cutadapt'] = "Primer removal process enabled"
if (params.figaro_enable) summary['FIGARO'] = "Optimizing microbiome rRNA gene trimming parameters for DADA2 enabled"
if (params.dbotu3_enable) summary['dbOTU3'] = "ASV clustering based on phylogeny, distribution and abundance enabled"
summary['Taxonomic database used'] = params.database
if (params.filter_table_by_tax_enable) summary['Filtering process'] = "Based on taxonomy enabled"
if (params.filter_table_by_data_enable) {
    if (params.filter_by_id && params.filter_by_frequency) {
        summary['Filtering processes'] = "Based on sample ID and frequency enabled"
    } else if (params.filter_by_id) {
        summary['Filtering process'] = "Based on sample ID enabled"
    } else {
        summary['Filtering process'] = "Based on frequency enabled"
    }
}
if (params.filter_contaminants_enable) summary['microDecon'] = "Sample decontamination process enabled"

log.info summary.collect { k,v -> "${k.padRight(24)}: $v" }.join("\n")
log.info "\033[1;34m-------------------------------------------------------------------\033[0m"

// Check the hostnames against configured profiles
checkHostname()

/*
 * VERIFY WORKFLOW VARIABLES
 */

/* Illumina workflow */
if (workflow.profile.contains('illumina')) {

    /* Verify Cutadapt parameters */
    if (params.cutadapt_enable) {
        if(params.primerF.isEmpty() || params.primerR.isEmpty()) {
            log.error "ERROR: no primer sequences have been provided. Please check and configure the '--primerF' and '--primerR' parameters in the illumina.config file"
            exit 1
        }
    }

    /* Verify FIGARO parameters */
    if (params.figaro_enable) {
        if(params.raw_data_dir.isEmpty()) {
            log.error "ERROR: raw data directory has not been configured. Please check and configure the '--raw_data_dir' parameter in the illumina.config file"
            exit 1
        }
        if(params.amplicon_length.isEmpty()) {
            log.error "ERROR: no expected amplicon size has been provided. Please check and configure the '--amplicon_length' parameter in the illumina.config file"
            exit 1
        }
    }

    /* Verify DADA2 parameters */
    if(params.FtrimLeft.isEmpty() || params.RtrimLeft.isEmpty() || params.FtruncLen.isEmpty() || params.RtruncLen.isEmpty() || params.truncQ.isEmpty() || params.FmaxEE.isEmpty() || params.RmaxEE.isEmpty() || params.pooling_method.isEmpty() || params.chimeras_method.isEmpty() ) {
        log.error "ERROR: DADA2 parameters have not been configured correctly. At least one of the parameters is not filled in. Please check and configure all parameters in the 'DADA2 process parameters' section of the illumina.config file"
        exit 1
    }

    /* Verify the taxonomic database */
    if(params.database.isEmpty()) {
        log.error "ERROR: No taxonomic database has been provided. Please check and configure the '--database' parameter in the illumina.config file"
        exit 1
    }

    /* Verify the list of taxa to exclude if params.filter_table_by_tax_enable process is activated */
    if (params.filter_table_by_tax_enable) {
        if (params.tax_to_exclude.isEmpty()) {
            log.error "ERROR: The list of taxa to exclude is empty. Please check and configure the '--tax_to_exclude' parameter in the illumina.config file"
            exit 1
        }
    }
    
    /* Verify parameters for filter_table_by_data process process if it is activated */
    if (params.filter_table_by_data_enable) {
        if (params.filter_by_id) {
            if (params.list_sample_to_remove.isEmpty()) {
                log.error "ERROR: The list of sample to remove is empty. Please check and configure the '--list_sample_to_remove' parameter in the illumina.config file"
                exit 1
            }
        } else {
            params.list_sample_to_remove = "none"
        }
        if (params.filter_by_frequency) {
            if (params.min_frequency_sample.isEmpty() || params.min_frequency_asv.isEmpty() || params.contingency_asv.isEmpty()) {
                log.error "ERROR: At least one of the frequency-based filter parameters is empty. Please check and configure the '--min_frequency_sample', '--min_frequency_asv' and/or '--contingency_asv' parameters in the illumina.config file"
                exit 1
            }
        } else {
            params.min_frequency_sample = "1"
            params.min_frequency_asv = "1"
            params.contingency_asv = "1"
        }
    }

    /* Verify parameters for filter_contaminants process process if it is activated */
    if (params.filter_contaminants_enable) {
        if (params.list_control_samples.isEmpty()) {
                log.error "ERROR: The list of control samples is empty. Please check and configure the '--control_samples' parameter in the illumina.config file"
                exit 1
        }
    }
}

/*
 *  SET UP WORKFLOW CHANNELS
 */

if (!workflow.profile.contains('test')) {
    channel
        .fromPath( params.excel_sample_file )
        .ifEmpty { error "ERROR: Cannot find the Sample Input Excel File at this path: ${params.excel_sample_file}. Please check and correct the parameter 'excel_sample_file' provided in the your analysis config file" }
        .set { sample_file }
}

/*
 * IMPORTING MODULES
 */

include { get_test_data } from './modules/get_test_data.nf'
include { excel2tsv } from './modules/excel2tsv.nf'
include { addpath_testdata } from './modules/excel2tsv.nf'
include { data_integrity } from './modules/data_integrity.nf'
include { q2_import_data } from './modules/qiime2.nf'
include { q2_cutadapt } from './modules/qiime2.nf'
include { figaro } from './modules/figaro.nf'
include { q2_dada2 } from './modules/qiime2.nf'
include { q2_dbOTU3 } from './modules/qiime2.nf'
include { q2_assign_taxo } from './modules/qiime2.nf'
include { q2_filter_table_by_tax } from './modules/qiime2.nf'
include { q2_filter_table_by_data } from './modules/qiime2.nf'
include { filter_contaminants } from './modules/filter_contaminants.nf'

/*
 * RUN MAIN WORKFLOW
 */

workflow {

    /*---------------------------------------*/
    /*                                       */
    /* General processes                     */
    /*                                       */
    /*---------------------------------------*/

    /* IF TEST WORKFLOW : Download raw data */
        if (workflow.profile.contains('test')) {
            get_test_data()
            ready = get_test_data.out.data4test_ready
            sample_file = get_test_data.out.xls
        } else {
            ready = channel.value('custom_workflow')
        }

    /* Convert Excel file */
        excel2tsv(sample_file, ready)
        if (workflow.profile.contains('test')) {
            addpath_testdata(excel2tsv.out.manifest_xls)
        }
        manifest = workflow.profile.contains('test') ? addpath_testdata.out.manifest : excel2tsv.out.manifest_xls

    /*---------------------------------------*/
    /*                                       */
    /* Illumina short reads analysis         */
    /*                                       */
    /*---------------------------------------*/

    if (!params.longreads && !params.stats_only && !params.dada2merge) {
    
        /* Verify data integrity */
            if (params.data_integrity_enable) {
                data_integrity(manifest,excel2tsv.out.metadata_xls)
            }

        /* Import data in QIIME2 format */
            q2_import_data(data_integrity.out.final_manifest)
        
        /* OPTIONAL: Primer removal using Cutadapt */
            if (params.cutadapt_enable) {
                q2_cutadapt(q2_import_data.out.imported_data)
            }
    
        /* OPTIONAL: Optimizing rRNA gene trimming parameters for DADA2 using FIGARO */
            if (params.figaro_enable) {
                figaro(ready)
            }

        /* ASV inference using DADA2 */
            dada2_input = params.cutadapt_enable ? q2_cutadapt.out.trimmed_data : q2_import_data.out.imported_data
            q2_dada2(dada2_input,excel2tsv.out.metadata_xls)

        /* ASV clustering using dbOTU3 */
            if (params.dbotu3_enable) {
                q2_dbOTU3(q2_dada2.out.dada2_table,q2_dada2.out.dada2_rep_seqs,excel2tsv.out.metadata_xls)
            }

        /* Taxonomic assignation of ASVs */
            asv_table = params.dbotu3_enable ? q2_dbOTU3.out.dbotu3_table : q2_dada2.out.dada2_table
            asv_sequences = params.dbotu3_enable ? q2_dbOTU3.out.dbotu3_seqs : q2_dada2.out.dada2_rep_seqs
            asv_outdir = params.dbotu3_enable ? q2_dbOTU3.out.dbotu3_outdir : q2_dada2.out.dada2_outdir
            asv_seqs_fasta = params.dbotu3_enable ? q2_dbOTU3.out.dbotu3_asv_seqs_fasta : q2_dada2.out.dada2_asv_seqs_fasta
            q2_assign_taxo(asv_sequences,asv_outdir)

        /* OPTIONAL: Filter ASV table and ASV sequences based on taxonomy */
            if (params.filter_table_by_tax_enable) {
                q2_filter_table_by_tax(asv_table,asv_sequences,q2_assign_taxo.out.taxonomy_assigned,q2_assign_taxo.out.taxonomy_tsv,excel2tsv.out.metadata_xls)
            }

        /* OPTIONAL: Filter ASV table and ASV sequences based on data */
            asv_table = params.filter_table_by_tax_enable ? q2_filter_table_by_tax.out.asv_table_tax_filtered_qza : asv_table
            asv_sequences = params.filter_table_by_tax_enable ? q2_filter_table_by_tax.out.asv_seqs_tax_filtered_qza : asv_sequences
            if (params.filter_table_by_data_enable) {
                q2_filter_table_by_data(asv_table,asv_sequences,excel2tsv.out.metadata_xls,q2_assign_taxo.out.taxonomy_tsv)
            }

        /* OPTIONAL: Filter contaminant ASVs using control samples */
            asv_table_tsv_contaminated = params.filter_table_by_data_enable ? q2_filter_table_by_data.out.final_asv_table_filtered_tsv : params.filter_table_by_tax_enable ? q2_filter_table_by_tax.out.asv_table_tax_filtered_tsv : q2_assign_taxo.out.asv_tax_table_tsv 
            asv_sequences_fasta_contaminated = params.filter_table_by_data_enable ? q2_filter_table_by_data.out.filter_table_by_data_seqs_fasta : params.filter_table_by_tax_enable ? q2_filter_table_by_tax.out.filter_table_by_tax_seqs_fasta : asv_seqs_fasta
            if (params.filter_contaminants_enable) {
                filter_contaminants(asv_table_tsv_contaminated,asv_sequences_fasta_contaminated,excel2tsv.out.metadata_xls)
            }

    }

}

/*
 * Completion notification
 */

workflow.onComplete {
    c_blue = "\033[1;34m";
    c_yellow = "\033[1;33m";
    c_green = "\033[1;32m";
    c_red = "\033[1;31m";
    c_reset = "\033[0m";

    if (workflow.success) {
        log.info """${c_blue}========================================================================${c_reset}
${c_yellow}SAMBA workflow${c_reset}: ${c_green}Pipeline completed successfully${c_reset}"""
    } else {
        checkHostname()
        log.info """${c_blue}========================================================================${c_reset}
${c_yellow}SAMBA workflow${c_reset}: ${c_red}Pipeline completed with errors${c_reset}"""
    }
}

/*
 * Other functions
 */

def SeBiMERHeader() {
    // Log colors ANSI codes
    c_red = '\033[1;31m'
    c_blue = '\033[1;34m'
    c_reset = '\033[0m'
    c_yellow = '\033[0;33m'
    c_purple = '\033[0;35m'

    return """    ${c_blue}------------------------------------------------------------------${c_reset}
    ${c_red}    __  __  __  .       __  __  ${c_reset}
    ${c_red}   \\   |_  |__) | |\\/| |_  |__)  ${c_reset}
    ${c_red}  __\\  |__ |__) | |  | |__ |  \\  ${c_reset}
                                            ${c_reset}
    ${c_yellow}  SAMBA workflow (${workflow.manifest.version})${c_reset}
                                            ${c_reset}
    ${c_purple}  Homepage: ${workflow.manifest.homePage}${c_reset}
    ${c_blue}-------------------------------------------------------------------${c_reset}
    """.stripIndent()
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;31m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;33m"
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
