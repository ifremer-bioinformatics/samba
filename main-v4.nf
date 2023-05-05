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
	--data_type			[str]	Set the type of your data. Can be: illumina, nanopore or pacbio.
	--singleEnd			[bool]	Set to true to specify that the inputs are single-end reads (default = false).
	--raw_read_length		[int]	Set the length of your raw reads e.g 250 or 300bp (default = 250).

	General parameters:
	--stat_var			[str]	Variable(s) of interest for statistical analyses (comma-separated list).
	--taxa_nb			[int]	Number of taxa to represent (default = "10").

	Other options:
	--outdir			[path]	The output directory where the results will be saved.
	-name 				[str]	Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
	--projectName			[str]	Name of the project.

	@@ PROCESS OPTIONS FOR ILLUMINA ANALYSIS @@

	Data integrity:
	--data_integrity_enable		[bool]	Data integrity checking step. Set to false to deactivate this step (default = true).
	--primer_filter			[int]	Percentage of primers supposed to be found in raw reads (default = 70).

	Cutadapt - primer removal:
	--cutadapt_enable		[bool]	Primer removal process. Set to false to deactivate this step (default = true).
	--primerF			[str]	Forward primer (to be used in Cutadapt cleaning step).
	--primerR			[str]	Reverse primer (to be used in Cutadapt cleaning step).
	--errorRate			[float]	Cutadapt error rate allowed to match primers (default = 0.1).

	FIGARO - optimizing trimming parameters for DADA2:
	--figaro_enable			[bool]	Process to identify optimal trimming parameters for DADA2. Set to false to deactivate this step (default = true).
	--raw_data_dir			[path]	Path to raw data directory.
	--amplicon_length		[str]	Length of expected amplicons.

	DADA2 - ASV inference:
	--FtrimLeft			[int]	The number of nucleotides to remove from the start of each forward read (default : 0 = no trimming).
	--RtrimLeft			[int]	The number of nucleotides to remove from the start of each reverse read (default : 0 = no trimming).
	--FtruncLen			[int]	Truncate forward reads after FtruncLen bases. Reads shorter than this are discarded (default : 0 = no trimming).
	--RtruncLen			[int]	Truncate reverse reads after RtruncLen bases. Reads shorter than this are discarded (default : 0 = no trimming).
	--FmaxEE			[float]	Forward reads with higher than max "expected errors" will be discarded (default = 2).
	--RmaxEE			[float]	Reverse with higher than max "expected errors" will be discarded (default = 2).
	--truncQ			[int]	Truncate reads at the first instance of a quality score less than or equal to minQ (default = 2).
	--n_read_learn			[int]	number of reads to use when training the error model (default = "1000000").
	--pooling_method		[str]	Method used to pool samples for denoising. Default = "independant". Set to "pseudo" if you want to approximate pooling of samples (see DADA2 documentation).
	--chimeras_method		[str]	Chimera detection method : default = "consensus". Set to "pooled" if the samples in the sequence table are all pooled together for bimera identification (see DADA2 documentation).

	swarm clustering of ASVs:
	--swarm_clustering_enable	[bool]	Set to true to activate the swarm ASV clustering process (default = false).

	dbOTU3 - Distribution based-clustering:
	--dbotu3_enable			[bool]	Distribution based-clustering step. Set to false to deactivate
 this step (default = true).
	--gen_crit			[float]	dbOTU3 Genetic criterion (default = 0.1).
	--abund_crit			[int]	dbOTU3 Abundance criterion (default = 10).
	--pval_crit			[float]	dbOTU3 P-value criterion (default = 0.0005).

	Taxonomic assignation:
	--database			[file]	Path to a trained Naive Bayes QIIME 2 classifier.
	--confidence			[float]	Confidence threshold for limiting taxonomic depth. Set to "disable" to disable confidence calculation, or 0 to calculate confidence but not apply it to limit the taxonomic depth of the assignments (default = 0.7).

	Filter ASV table and ASV sequences based on taxonomy:
	--filter_table_by_tax_enable	[bool]	Set to true to filter ASV table and ASV sequences based on taxonomic assignation (default = false).
	--filtering_type		[str]	Type of filtering: 'exclude' if you want to exclude some taxa or 'include' if you want to keep only some taxa
	--tax_to_filter			[str]	List of taxa you want to filter (comma-separated list).

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

	Differential abundance testing using ANCOM-BC:
	--ancombc_enable		[bool]	Set to true to run differential abundance analysis on your data (default = false).
	--use_custom_reference		[bool]	Set to true if you have a reference group for the differential abundance analysis (default = false).
	--reference_level		[str]	Group of samples to use as reference. Syntax: "column_name::column_value"
	--ancombc_formula		[str]	Variable to analyse for differential abundance. Can be: "column_name" (test only for difference between values in this column) ; "column_name_1,column_name2,..." (as previously but for multiple variable) or "column_name_1+column_name2" (test difference using interacting variables).
	--p_adj_method			[str]	Pvalue adjusting method. Can be: 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none' (default = "holm").
	--max_iter			[int]	Maximum number of iterations for the E-M algorithm (default = "100").
	--alpha				[float]	Level of significance (default = "0.05").

	Functional predictions using PICRUSt2:
	--picrust2_enable		[bool]  Set to true to run functional predictions (only efficient for the 16S marker gene) (default = false).
	--traits_db			[str]	Gene families to predict. Can be: COG, EC, KO, PFAM, TIGRFAM (default = "EC,KO").
	--nsti				[float]	Max nsti value accepted. NSTI cut-off of 2 should eliminate junk sequences (default = "0.2").
	--hsp_method			[str]	HSP method of your choice. 'mp' the most accurate prediction method, faster method: 'pic' (default = "mp").
	--min_reads			[int]	Minimum number of reads across all samples for each input ASV (default = "1").
	--min_samples			[int]	Minimum number of samples that an ASV needs to be identfied within (default = "1").
	--picrust2_tested_variable	[str]	Variable(s) of interest for functional predictions (comma-separated list).

	@@ PROCESS OPTIONS FOR NANOPORE ANALYSIS @@

	Length filtering of raw reads:
	--nanopore_read_minlength	[int]	Minimum length of raw nanopore reads.
	--nanopore_read_maxlength	[int]	Maximum length of raw nanopore reads.

	Read mapping using minimap2:
	--batch_size			[int]	Minibatch size for mapping (default = "25M").
	--minimap2_preset		[str]	minimap2 preset (default = "map-ont").
	--minimap2_db			[path]	Path to a minimap2 index format of the taxonomic database.

	Get count table:
	--ref_tax			[path]	Path to the tabulated file containing the taxonomic reference ID and the associated taxonomy
	--tax_rank			[int]	Taxonomic level to analyse data. Can be: Kingdom;Phylum;Class;Order;Genus;Species.
	--kingdom			[str]	Kingdom studied 

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
if (workflow.profile == 'illumina,singularity') {
    illumina_params_file = file("${baseDir}/conf/illumina.config", checkIfExists: true)
    illumina_params_file.copyTo("${params.outdir}/00_pipeline_config/illumina.config")
}
if (workflow.profile == 'illumina_test,singularity') {
    illumina_test_params_file = file("${baseDir}/conf/illumina_test.config", checkIfExists: true)
    illumina_test_params_file.copyTo("${params.outdir}/00_pipeline_config/illumina_test.config")
}
if (workflow.profile == 'nanopore_test,singularity') {
   nanopore_test_params_file = file("${baseDir}/conf/nanopore_test.config", checkIfExists: true)
   nanopore_test_params_file.copyTo("${params.outdir}/00_pipeline_config/nanopore_test.config")
}
if (workflow.profile == 'nanopore,singularity') {
   nanopore_params_file = file("${baseDir}/conf/nanopore.config", checkIfExists: true)
   nanopore_params_file.copyTo("${params.outdir}/00_pipeline_config/nanopore.config")
}

file_excel_sample_file = new File(params.excel_sample_file)
if (params.data_type == 'illumina') file_database = new File(params.database)
if (params.data_type == 'nanopore') file_nanopore_database = new File(params.minimap2_db)
if (params.data_type == 'nanopore') file_nanopore_ref_tax = new File(params.ref_tax)

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
if (workflow.profile == 'illumina,singularity') summary['Data type'] = 'Illumina short reads'
if (workflow.profile == 'illumina_test,singularity') summary['Data type'] = 'Illumina short reads test workflow'
if (workflow.profile == 'nanopore,singularity') summary['Data type'] = 'Nanopore long reads'
if (workflow.profile == 'nanopore_test,singularity') summary['Data type'] = 'Nanopore long reads test workflow'
summary['Sample Input Excel File'] = file_excel_sample_file.name
if (params.data_type == 'illumina') {
    summary['Data integrity'] = params.data_integrity_enable ? "Enabled" : "Disabled"
    summary['Primer removal using Cutadapt'] = params.cutadapt_enable ? "Enabled" : "Disabled"
    summary['Optimizing DADA2 parameters using FIGARO'] = params.figaro_enable ? "Enabled" : "Disabled"
    summary['ASV clustering using swarm'] = params.swarm_clustering_enable ? "Enabled" : "Disabled"
    summary['ASV clustering using dbOTU3'] = params.dbotu3_enable ? "Enabled" : "Disabled"
    summary['Taxonomic database used'] = file_database.name
    summary['Taxonomy filtering'] = params.filter_table_by_tax_enable ? "Enabled" : "Disabled"
    if (params.filter_table_by_tax_enable) {
        summary['    |_ Type of tax filtering'] = params.filtering_type
        summary['    |_ Taxa to filter'] = params.tax_to_filter
    }
    summary['Data filtering'] = params.filter_table_by_data_enable ? "Enabled" : "Disabled"
    if (params.filter_table_by_data_enable) {
        if (params.filter_by_id && params.filter_by_frequency) {
            summary['    |_ based on sample ID and frequency'] = "Enabled"
        } else if (params.filter_by_id) {
            summary['    |_ based on sample ID'] = "Enabled"
        } else {
            summary['    |_ based on frequency'] = "Enabled"
        }
    }
    summary['Sample decontamination (microDecon)'] = params.filter_contaminants_enable ? "Enabled" : "Disabled"
    summary['Differential abundance testing (ANCOM-BC)'] = params.ancombc_enable ? "Enabled" : "Disabled"
    summary['Functional predictions (PICRUSt2)'] = params.picrust2_enable ? "Enabled" : "Disabled"
    if (params.picrust2_enable) {
        summary['    |_ Predicted Gene Families'] = params.traits_db
    }
}
if (params.data_type == 'nanopore') {
    summary['Taxonomic database used'] = file_nanopore_database.name
    summary['Reference taxonomy file'] = file_nanopore_ref_tax.name
    summary['Taxonomic level inspected'] = params.tax_rank
}

log.info summary.collect { k,v -> "${k.padRight(42)}: $v" }.join("\n")
log.info "\033[1;34m-------------------------------------------------------------------\033[0m"

// Check the hostnames against configured profiles
checkHostname()

/*
 * VERIFY WORKFLOW VARIABLES
 */

/* Illumina workflow */
if (params.data_type == 'illumina') {

    /* Verify Sample file */
    if (!workflow.profile.contains('test')) {
        if(params.excel_sample_file.isEmpty()) {
            log.error "ERROR: Cannot find the Sample Input Excel File at this path: ${params.excel_sample_file}. Please check and correct the parameter 'excel_sample_file' provided in the your analysis config file"
            exit 1
        }
    }

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
    if(params.FtrimLeft.isEmpty() || params.RtrimLeft.isEmpty() || params.FtruncLen.isEmpty() || params.RtruncLen.isEmpty() || params.truncQ.isEmpty() || params.FmaxEE.isEmpty() || params.RmaxEE.isEmpty() || params.n_read_learn.isEmpty() || params.pooling_method.isEmpty() || params.chimeras_method.isEmpty() ) {
        log.error "ERROR: DADA2 parameters have not been configured correctly. At least one of the parameters is not filled in. Please check and configure all parameters in the 'DADA2 process parameters' section of the illumina.config file"
        exit 1
    }

    /* Set dbOTU3 to false if swarm clustering is activated */
    if (params.swarm_clustering_enable) {
        params.dbotu3_enable = false
    }

    /* Verify the taxonomic database */
    if(params.database.isEmpty()) {
        log.error "ERROR: No taxonomic database has been provided. Please check and configure the '--database' parameter in the illumina.config file"
        exit 1
    }

    /* Verify the list of taxa to filter if params.filter_table_by_tax_enable process is activated */
    if (params.filter_table_by_tax_enable) {
        if (params.filtering_type.isEmpty() || params.tax_to_filter.isEmpty()) {
            log.error "ERROR: Type of tax filtering and/or the list of taxa to filter is empty. Please check and configure the '--filtering_type' and '--tax_to_filter' parameters in the illumina.config file"
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

    /* Verify parameters for filter_contaminants process if it is activated */
    if (params.filter_contaminants_enable) {
        if (params.list_control_samples.isEmpty()) {
                log.error "ERROR: The list of control samples is empty. Please check and configure the '--control_samples' parameter in the illumina.config file"
                exit 1
        }
    }

    /* Verify parameters for ANCOM-BC if it is activated */
    if (params.ancombc_enable) {
        if (params.ancombc_formula.isEmpty() || params.p_adj_method.isEmpty() || params.max_iter.isEmpty() || params.alpha.isEmpty()) {
            log.error "ERROR: At least one of the ANCOM-BC parameters is empty. Please check and configure the '--ancombc_formula', '--p_adj_method', '--max_iter' and/or '--alpha' parameters in the illumina.config file"
            exit 1
        }
        if (params.use_custom_reference) {
            if (params.reference_level.isEmpty()) {
                log.error "ERROR: No condition is provided to used as reference for the differential abundance analysis (ANCOM-BC). Please check and configure the '--reference_level' parameter in the illumina.config file"
                exit 1
            }
        }
    }

    /* Verify parameters for functional predictions if it is activated */
    if (params.picrust2_enable) {
        if (params.nsti.isEmpty() || params.hsp_method.isEmpty() || params.min_reads.isEmpty() || params.min_samples.isEmpty() || params.picrust2_tested_variable.isEmpty()) {
            log.error "ERROR: At least one of the parameters for the functional predictions process is empty. Please check and configure the '--nsti', '--hsp_method', '--min_reads', '--min_samples' and/or '--picrust2_tested_variable' parameters in the illumina.config file"
            exit 1
        }
    }
}

/* Nanopore workflow */
if (params.data_type == 'nanopore') {
    if (params.nanopore_read_minlength.isEmpty() || params.nanopore_read_maxlength.isEmpty()) {
        log.error "ERROR: No miminum and/or maximum length for raw nanopore reads has/have been provided. Please check and configure the '--nanopore_read_minlength' and/or '--nanopore_read_maxlength' parameters in the nanopore.config file"
        exit 1
    }
    if (params.batch_size.isEmpty() || params.minimap2_preset.isEmpty() || params.minimap2_db.isEmpty()) {
        log.error "ERROR: Please check and configure all the minimap2 parameters in the nanopore.config file: '--batch_size', '--minimap2_preset' and '--minimap2_db'"
        exit 1
    }
    if (params.ref_tax.isEmpty() || params.tax_rank.isEmpty()) {
        log.error "ERROR: No reference taxonomy file and/or taxonomic level to analyse has/have been provided. Please check and configure the '--ref_tax' and/or '--tax_rank' parameters in the nanopore.config file"
        exit 1
    }
    if (params.kingdom.isEmpty()) {
        log.error "ERROR: No Kingdom studied has been provided. Please check and configure the '--kingdom' parameter in the nanopore.config file"
        exit 1
    }
}

/*
 *  SET UP WORKFLOW CHANNELS
 */

if (!workflow.profile.contains('test')) {
    channel
        .fromPath( params.excel_sample_file )
        .set { sample_file }
}

if (params.data_type == 'illumina') {
    if (params.ancombc_enable) {
        channel
            .from(params.ancombc_formula)
            .splitCsv(sep : ',', strip : true)
            .flatten()
            .set { ancombc_formula_ch }
    }
    
    if (params.picrust2_enable) {
        channel
            .from(params.picrust2_tested_variable)
            .splitCsv(sep : ',', strip : true)
            .flatten()
            .set { picrust2_tested_variable_ch }
    }
}

channel
    .from(params.stat_var)
    .splitCsv(sep : ',', strip : true)
    .flatten()
    .set { stat_var_ch }

/*
 * IMPORTING MODULES
 */

/* Illumina modules */
include { get_test_data } from './modules/get_test_data.nf'
include { excel2tsv } from './modules/excel2tsv.nf'
include { addpath_testdata } from './modules/excel2tsv.nf'
include { data_integrity } from './modules/data_integrity.nf'
include { q2_import_data } from './modules/qiime2.nf'
include { q2_cutadapt } from './modules/qiime2.nf'
include { figaro } from './modules/figaro.nf'
include { q2_dada2 } from './modules/qiime2.nf'
include { swarm_clustering_processing } from './modules/swarm_clustering.nf'
include { swarm_clustering_format_output } from './modules/swarm_clustering.nf'
include { q2_dbOTU3 } from './modules/qiime2.nf'
include { q2_assign_taxo } from './modules/qiime2.nf'
include { q2_filter_table_by_tax } from './modules/qiime2.nf'
include { q2_filter_table_by_data } from './modules/qiime2.nf'
include { filter_contaminants } from './modules/filter_contaminants.nf'
include { q2_asv_phylogeny } from './modules/qiime2.nf'
include { format_final_outputs } from './modules/format_final_outputs.nf'
include { q2_ancombc } from './modules/qiime2.nf'
include { picrust2 } from './modules/picrust2.nf'
include { create_phyloseq } from './modules/R.nf'

/* Nanopore modules */
include { nanopore_read_length_filter } from './modules/nanopore.nf'
include { nanopore_mapping } from './modules/nanopore.nf'
include { nanopore_getfasta } from './modules/nanopore.nf'
include { nanopore_count_table } from './modules/nanopore.nf'
include { nanopore_phyloseq_obj } from './modules/R.nf'
include { nanopore_alpha_diversity } from './modules/R.nf'

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

    if (params.data_type == "illumina") {
    
        /* Verify data integrity */
            if (params.data_integrity_enable) {
                data_integrity(manifest,excel2tsv.out.metadata_xls)
            }
        
        /* Import data in QIIME2 format */
            q2_import_data(manifest)
        
        /* OPTIONAL: Primer removal using Cutadapt */
            if (params.cutadapt_enable) {
                q2_cutadapt(q2_import_data.out.imported_data)
            }
    
        /* OPTIONAL: Optimizing rRNA gene trimming parameters for DADA2 using FIGARO */
            if (params.figaro_enable) {
                figaro(ready)
                optimal_values = figaro.out.figaro_csv.splitCsv(header: true, sep:',', limit:1).map { row -> tuple( row."trimPosition_R1", row."trimPosition_R2", row."maxExpectedError_R1", row."maxExpectedError_R2") }
            }

        /* ASV inference using DADA2 */
            /* ~~~ input management ~~~ */
            dada2_input = params.cutadapt_enable ? q2_cutadapt.out.trimmed_data : q2_import_data.out.imported_data
            dada2_trimParams = params.figaro_enable ? optimal_values : Channel.from(params.FtruncLen,params.RtruncLen,params.FmaxEE,params.RmaxEE).collect()
            /* ~~~ process ~~~ */
            q2_dada2(dada2_input,excel2tsv.out.metadata_xls,dada2_trimParams)

        /* OPTIONAL: swarm clustering */
            if (params.swarm_clustering_enable) {
                swarm_clustering_processing(q2_dada2.out.dada2_asv_seqs_fasta_abundance)
                swarm_clustering_format_output(q2_dada2.out.dada2_table_tsv,swarm_clustering_processing.out.asv_swarm_cluster_list_tsv,swarm_clustering_processing.out.asv_swarm_seqs_clusters_fasta,excel2tsv.out.metadata_xls)
            }

        /* ASV clustering using dbOTU3 */
            if (params.dbotu3_enable && !params.swarm_clustering_enable) {
                q2_dbOTU3(q2_dada2.out.dada2_table,q2_dada2.out.dada2_rep_seqs,excel2tsv.out.metadata_xls)
            }

        /* Taxonomic assignation of ASVs */
            /* ~~~ input management ~~~ */
            asv_table = params.swarm_clustering_enable ? swarm_clustering_format_output.out.swarm_asv_table : params.dbotu3_enable ? q2_dbOTU3.out.dbotu3_table : q2_dada2.out.dada2_table
            asv_sequences = params.swarm_clustering_enable ? swarm_clustering_format_output.out.swarm_asv_seqs : params.dbotu3_enable ? q2_dbOTU3.out.dbotu3_seqs : q2_dada2.out.dada2_rep_seqs
            asv_outdir = params.swarm_clustering_enable ? swarm_clustering_format_output.out.swarm_asv_outdir : params.dbotu3_enable ? q2_dbOTU3.out.dbotu3_outdir : q2_dada2.out.dada2_outdir
            asv_sequences_fasta = params.swarm_clustering_enable ? swarm_clustering_processing.out.asv_swarm_seqs_clusters_fasta : params.dbotu3_enable ? q2_dbOTU3.out.dbotu3_asv_seqs_fasta : q2_dada2.out.dada2_asv_seqs_fasta
            /* ~~~ process ~~~ */
            q2_assign_taxo(asv_sequences,asv_outdir)

        /* OPTIONAL: Filter contaminant ASVs using control samples */
            if (params.filter_contaminants_enable) {
                filter_contaminants(q2_assign_taxo.out.asv_tax_table_tsv,asv_sequences_fasta,excel2tsv.out.metadata_xls)
            }

        /* OPTIONAL: Filter ASV table and ASV sequences based on taxonomy */
            if (params.filter_table_by_tax_enable) {
                /* ~~~ input management ~~~ */
                asv_table = params.filter_contaminants_enable ? filter_contaminants.out.decontam_ASV_table_qza : asv_table
                asv_sequences = params.filter_contaminants_enable ? filter_contaminants.out.decontam_ASV_seqs_qza : asv_sequences
                /* ~~~ process ~~~ */
                q2_filter_table_by_tax(asv_table,asv_sequences,q2_assign_taxo.out.taxonomy_assigned,q2_assign_taxo.out.taxonomy_tsv,excel2tsv.out.metadata_xls)
            }

        /* OPTIONAL: Filter ASV table and ASV sequences based on data */
            if (params.filter_table_by_data_enable) {
                /* ~~~ input management ~~~ */
                asv_table = params.filter_table_by_tax_enable ? q2_filter_table_by_tax.out.asv_table_tax_filtered_qza : asv_table
                asv_sequences = params.filter_table_by_tax_enable ? q2_filter_table_by_tax.out.asv_seqs_tax_filtered_qza : asv_sequences
                /* ~~~ process ~~~ */
                q2_filter_table_by_data(asv_table,asv_sequences,excel2tsv.out.metadata_xls,q2_assign_taxo.out.taxonomy_tsv)
            }

        /* Construct a phylogeny of ASVs */
            /* ~~~ input management ~~~ */
            asv_seqs_phylo = params.filter_table_by_data_enable ? q2_filter_table_by_data.out.final_asv_seqs_filtered_qza : asv_sequences
            /* ~~~ process ~~~ */
            q2_asv_phylogeny(asv_seqs_phylo) 

        /* Format final outputs */
            /* ~~~ input management ~~~ */
            final_asv_table_tsv = params.filter_table_by_data_enable ? q2_filter_table_by_data.out.final_asv_table_filtered_tsv : params.filter_table_by_tax_enable ? q2_filter_table_by_tax.out.asv_table_tax_filtered_tsv : params.filter_contaminants_enable ? filter_contaminants.out.decontam_ASV_table_tsv : q2_assign_taxo.out.asv_tax_table_tsv
            final_asv_sequences_fasta = params.filter_table_by_data_enable ? q2_filter_table_by_data.out.filter_table_by_data_seqs_fasta : params.filter_table_by_tax_enable ? q2_filter_table_by_tax.out.filter_table_by_tax_seqs_fasta : params.filter_contaminants_enable ? filter_contaminants.out.decontam_ASV_seqs_fasta : asv_sequences_fasta
            final_asv_table_biom = params.filter_table_by_data_enable ? q2_filter_table_by_data.out.final_asv_table_filtered_biom : params.filter_table_by_tax_enable ? q2_filter_table_by_tax.out.asv_table_tax_filtered_biom : params.filter_contaminants_enable ? filter_contaminants.out.decontam_ASV_table_biom : q2_assign_taxo.out.asv_tax_table_biom
            final_asv_table_qza = params.filter_table_by_data_enable ? q2_filter_table_by_data.out.final_asv_table_filtered_qza : params.filter_table_by_tax_enable ? q2_filter_table_by_tax.out.asv_table_tax_filtered_qza : params.filter_contaminants_enable ? filter_contaminants.out.decontam_ASV_table_qza : asv_table
            final_asv_seqs_qza = params.filter_table_by_data_enable ? q2_filter_table_by_data.out.final_asv_seqs_filtered_qza : params.filter_table_by_tax_enable ? q2_filter_table_by_tax.out.asv_seqs_tax_filtered_qza : params.filter_contaminants_enable ? filter_contaminants.out.decontam_ASV_seqs_qza : asv_sequences
            /* ~~~ process ~~~ */
            format_final_outputs(final_asv_table_tsv,final_asv_sequences_fasta,final_asv_table_biom,final_asv_table_qza,final_asv_seqs_qza)

        /* OPTIONAL: Differential abundance testing using ANCOM-BC */
            if (params.ancombc_enable) {
                q2_ancombc(final_asv_table_qza,excel2tsv.out.metadata_xls,q2_assign_taxo.out.taxonomy_assigned,ancombc_formula_ch)
            }

        /* OPTIONAL: Functional predictions using PICRUSt2 */
        if (params.picrust2_enable) {
            picrust2(final_asv_table_biom,final_asv_sequences_fasta)
        }

        /* Create the phyloseq object for statistical analyses */
            create_phyloseq(final_asv_table_tsv,excel2tsv.out.metadata_xls,q2_asv_phylogeny.out.asv_phylogeny_nwk)

    }

    /*---------------------------------------*/
    /*                                       */
    /* Nanopore long reads analysis          */
    /*                                       */
    /*---------------------------------------*/

    if (params.data_type == 'nanopore') {

        /* Raw reads length filter */
            /* ~~~ input management ~~~ */
            nanopore_manifest = manifest.splitCsv(header: true, sep:'\t').map { row -> tuple( row."sample-id", file(row."absolute-filepath")) }
            /* ~~~ process ~~~ */
            nanopore_read_length_filter(nanopore_manifest)

        /* Mapping of nanopore reads against tax database using minimap2 */
            nanopore_mapping(nanopore_read_length_filter.out.filtered_nanopore_fastq)

        /* Collect all Nanopore reads to a FASTA file */
            nanopore_getfasta(nanopore_read_length_filter.out.filtered_nanopore_fastq)
            nanopore_fasta = nanopore_getfasta.out.nanopore_sequences_fasta.collectFile(name : 'nanopore_sequences_all_samples.fasta', newLine : false, storeDir : "${params.outdir}/${params.nanopore_getfasta_results}").subscribe { println "All Nanopore sequences contained in the sample set are saved to the FASTA file : $it" }

        /* Get the count table */
            nanopore_count_table(nanopore_mapping.out.nanopore_mapped_reads.collect())

        /* Create the phyloseq object for statistical analyses */
           nanopore_phyloseq_obj(nanopore_count_table.out.nanopore_count_table,excel2tsv.out.metadata_xls)

        /* Run alpha diversity analyses */
           nanopore_alpha_diversity(nanopore_phyloseq_obj.out.phy_obj.collect(),stat_var_ch)

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
