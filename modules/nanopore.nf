process nanopore_read_length_filter {

    tag "$sample"
    label 'nanopore_env'
    label 'multithreads'

    publishDir "${params.outdir}/${params.nanopore_read_length_filtering_results}", mode: 'copy', pattern: '*.filtered.fastq.gz'
    publishDir "${params.outdir}/${params.report_dirname}/98_version", mode: 'copy', pattern: 'v_seqkit.txt'

    input:
        tuple val(sample), path(nanopore_fastq)

    output:
        tuple val(sample), path('*.filtered.fastq.gz'), emit: filtered_nanopore_fastq

    script:
    """
    NANOPORE_01_read_length_filter.sh ${task.cpus} ${params.nanopore_read_minlength} ${params.nanopore_read_maxlength} ${nanopore_fastq} ${sample}_read_length_filter.cmd &> ${sample}_read_length_filter.log 2>&1
    seqkit version | cut -d ' ' -f2 > v_seqkit.txt
    """
}

process nanopore_mapping {

    tag "$sample"
    label 'nanopore_env'
    label 'medRAM'

    publishDir "${params.outdir}/${params.nanopore_mapping_results}", mode: 'copy', pattern: '*.bam'
    publishDir "${params.outdir}/${params.report_dirname}/98_version", mode: 'copy', pattern: 'v_minimap2.txt'

    input:
        tuple val(sample), path(filtered_nanopore_fastq)

    output:
       path('*.bam'), emit: nanopore_mapped_reads
       path('v_minimap2.txt')

    script:
    """
    NANOPORE_02_mapping.sh ${task.cpus} ${params.batch_size} ${params.minimap2_preset} ${params.minimap2_db} ${filtered_nanopore_fastq} ${sample}.bam ${sample}_samtools_sort.log ${sample}_minimap2_mapping.cmd &> ${sample}_minimap2_mapping.log 2>&1
    minimap2 --version > v_minimap2.txt
    """
}

process nanopore_getfasta {

    tag "$sample"
    label 'qiime2_env'

    publishDir "${params.outdir}/${params.nanopore_getfasta_results}", mode: 'copy', pattern: '*_sequences.fasta'
    publishDir "${params.outdir}/${params.report_dirname}/98_version", mode: 'copy', pattern: 'v_seqtk.txt'

    input:
        tuple val(sample), path(fastq)
         
    output:
        path('*_sequences.fasta'), emit: nanopore_sequences_fasta

    script:
    """
    seqtk seq -a ${fastq} > ${sample}_sequences.fasta
    echo '1.3-r106' > v_seqtk.txt
    """

}

process nanopore_count_table {

    label 'biopython'

    publishDir "${params.outdir}/${params.nanopore_count_table_results}", mode: 'copy', pattern: '*.tsv'

    input:
        path(bam)

    output:
        path('samples.tsv'), emit: nanopore_count_table

    script:
    """
    NANOPORE_03_count_table.py -b . -t ${params.ref_tax} -r ${params.tax_rank} -o samples.tsv &> nanopore_count_table.log 2>&1
    """
}
