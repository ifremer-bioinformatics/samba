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

    The typical command for running the pipeline after filling the conf/custom.config file is as follows:

	nextflow run main.nf -profile <shortreadstest/longreadstest/custom>,singularity [-c <institute_config_file>]

	Mandatory:
	--excel_sample_file	[path]	Path to the XLS input file (EXCEL 97-2004) containing the manifest and metadata sheets.
	--longreads		[bool]	Set to true to specify that the inputs are long reads (Nanopore/Pacbio) (default = false for illumina short reads).
	--singleEnd		[bool]	Set to true to specify that the inputs are single-end reads (default = false).

	Other options:
	--outdir		[path]	The output directory where the results will be saved.
	-name 			[str]	Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
	--projectName		[str]	Name of the project.

	@@ PROCESS OPTIONS @@

	Data integrity:
	--data_integrity_enable	[bool]	Data integrity checking step. Set to false to deactivate this step (default = true).
	--primer_filter		[str]	Percentage of primers supposed to be found in raw reads (default = 70).

	Cutadapt - primer removal:
	--cutadapt_enable	[bool]	Primer removal process. Set to false to deactivate this step (default = true).
	--primerF		[str]	Forward primer (to be used in Cutadapt cleaning step).
	--primerR		[str]	Reverse primer (to be used in Cutadapt cleaning step).
	--errorRate		[str]	Cutadapt error rate allowed to match primers (default = 0.1).

	FIGARO - optimizing trimming parameters for DADA2:
	--figaro_enable		[bool]	Process to identify optimal trimming parameters for DADA2. Set to false to deactivate this step (default = true).
	--raw_data_dir		[path]	Path to raw data directory.
	--amplicon_length	[str]	Length of expected amplicons.

	DADA2 - ASV inference:
	--FtrimLeft		[str]	The number of nucleotides to remove from the start of each forward read (default : 0 = no trimming).
	--RtrimLeft		[str]	The number of nucleotides to remove from the start of each reverse read (default : 0 = no trimming).
	--FtruncLen		[str]	Truncate forward reads after FtruncLen bases. Reads shorter than this are discarded (default : 0 = no trimming).
	--RtruncLen		[str]	Truncate reverse reads after RtruncLen bases. Reads shorter than this are discarded (default : 0 = no trimming).
	--FmaxEE		[str]	Forward reads with higher than max "expected errors" will be discarded (default = 2).
	--RmaxEE		[str]	Reverse with higher than max "expected errors" will be discarded (default = 2).
	--truncQ		[str]	Truncate reads at the first instance of a quality score less than or equal to minQ (default = 2).
	--pooling_method	[str]	Method used to pool samples for denoising. Default = "independant". Set to "pseudo" if you want to approximate pooling of samples (see DADA2 documentation).
	--chimeras_method	[str]	Chimera detection method : default = "consensus". Set to "pooled" if the samples in the sequence table are all pooled together for bimera identification (see DADA2 documentation).

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
if (workflow.profile.contains('custom')) {
    custom_params_file = file("${baseDir}/conf/custom.config", checkIfExists: true)
    custom_params_file.copyTo("${params.outdir}/00_pipeline_config/custom.config")
}
if (workflow.profile.contains('shortreadstest')) {
    shortreadstest_params_file = file("${baseDir}/conf/shortreadstest.config", checkIfExists: true)
    shortreadstest_params_file.copyTo("${params.outdir}/00_pipeline_config/shortreadstest.config")
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
if (workflow.profile.contains('custom')) summary['Data type'] = params.longreads ? 'Long reads' : 'Illumina short reads'
if (workflow.profile.contains('test')) summary['Data type'] = params.longreads ? 'Long reads test workflow' : 'Illumina short reads test workflow'
summary['Sample Input Excel File'] = params.excel_sample_file
if (params.data_integrity_enable) summary['Data integrity'] = "Data integrity checking process enabled"
if (params.cutadapt_enable) summary['Cutadapt'] = "Primer removal process enabled"
if (params.figaro_enable) summary['FIGARO'] = "Optimizing microbiome rRNA gene trimming parameters for DADA2 enabled"

log.info summary.collect { k,v -> "${k.padRight(24)}: $v" }.join("\n")
log.info "\033[1;34m-------------------------------------------------------------------\033[0m"

// Check the hostnames against configured profiles
checkHostname()

/*
 * VERIFY WORKFLOW VARIABLES
 */

if (!workflow.profile.contains('custom')) {

    /* Verify Cutadapt parameters */
    if (params.cutadapt_enable) {
        if(params.primerF.isEmpty() || params.primerR.isEmpty()) {
            log.error "ERROR: no primer sequences have been provided. Please check and configure the '--primerF' and '--primerR' parameters in the custom.config file"
            exit 1
        }
    }

    /* Verify FIGARO parameters */
    if (params.figaro_enable) {
        if(params.raw_data_dir.isEmpty()) {
            log.error "ERROR: raw data directory has not been configured. Please check and configure the '--raw_data_dir' parameter in the custom.config file"
            exit 1
        }
        if(params.amplicon_length.isEmpty()) {
            log.error "ERROR: no expected amplicon size has been provided. Please check and configure the '--amplicon_length' parameter in the custom.config file"
            exit 1
        }
    }

    /* Verify DADA2 parameters */
    if(params.FtrimLeft.isEmpty() || params.RtrimLeft.isEmpty() || params.FtruncLen.isEmpty() || params.RtruncLen.isEmpty() || params.truncQ.isEmpty() || params.FmaxEE.isEmpty() || params.RmaxEE.isEmpty() || params.pooling_method.isEmpty() || params.chimeras_method.isEmpty() ) {
        log.error "ERROR: DADA2 parameters have not been configured correctly. At least one of the parameters is not filled in. Please check and configure all paramters in the 'DADA2 process parameters' section of the custom.config file"
        exit 1
    }
}

/*
 *  SET UP WORKFLOW CHANNELS
 */

if (!workflow.profile.contains('test')) {
    channel
        .fromPath( params.excel_sample_file )
        .ifEmpty { error "ERROR: Cannot find the Sample Input Excel File at this path: ${params.excel_sample_file}. Please check and correct the parameter 'excel_sample_file' provided in the custom.config file" }
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

    if (!params.longreads) {
    
        /* Verify data integrity */
            if (params.data_integrity_enable && !params.stats_only && !params.dada2merge) {
                data_integrity(manifest,excel2tsv.out.metadata_xls)
            }

        /* Import data in QIIME2 format */
            if (!params.stats_only && !params.dada2merge) {
                q2_import_data(data_integrity.out.final_manifest)
            }
        
        /* OPTIONAL: Primer removal using Cutadapt */
            if (!params.stats_only && !params.dada2merge) {
                q2_cutadapt(q2_import_data.out.imported_data)
            }
    
        /* OPTIONAL: Optimizing rRNA gene trimming parameters for DADA2 using FIGARO */
            if (!params.stats_only && !params.dada2merge) {
                figaro(ready)
            }

        /* ASV inference using DADA2 */
            dada2_input = params.cutadapt_enable ? q2_cutadapt.out.trimmed_data : q2_import_data.out.imported_data
            if (!params.stats_only && !params.dada2merge) {
                q2_dada2(dada2_input,excel2tsv.out.metadata_xls)
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
