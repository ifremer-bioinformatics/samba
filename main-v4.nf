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
	--outdir [path]			The output directory where the results will be saved.
	-name [str]			Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
	--projectName [str]		Name of the project.

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

log.info summary.collect { k,v -> "${k.padRight(24)}: $v" }.join("\n")
log.info "\033[1;34m-------------------------------------------------------------------\033[0m"

// Check the hostnames against configured profiles
checkHostname()

/*
 * VERIFY WORKFLOW VARIABLES
 */


/*
 *  SET UP WORKFLOW CHANNELS
 */

if (!workflow.profile.contains('test')) {
    channel
        .fromPath( params.excel_sample_file )
        .ifEmpty { error "Cannot find the Sample Input Excel File at this path: ${params.excel_sample_file}. Please check and correct the parameter 'excel_sample_file' provided in the custom.config file" }
        .set { sample_file }
}

/*
 * IMPORTING MODULES
 */

include { get_test_data } from './modules/get_test_data.nf'
include { excel2tsv } from './modules/excel2tsv.nf'
include { addpath_testdata } from './modules/excel2tsv.nf'
include { data_integrity } from './modules/data_integrity.nf'
include { import } from './modules/qiime2.nf'

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
                import(data_integrity.out.final_manifest)
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
