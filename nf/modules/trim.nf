/*
 * Modules for read trimming and other manipulations
 */

nextflow.enable.dsl=2

/*
 * This module defines a process for trimming a single-end fastQ file
 *
 * Config-defined parameters:
 * --------------------------
 *   + adapter_seq: the adapter sequence
 *   + trim_qual: the quality score for trimming
 *   + trim_ncores: the number of cores to use for trimming
 */
process trim_fq_single {
  conda params.HOME_REPO + '/nf/envs/cutadapt.yaml'
  input:
    tuple val(sequence_id), file(fastq_file)

  output:
    tuple val(sequence_id), file(trimmed_reads), file(trim_report)

  script:
    trimmed_reads = "${sequence_id}_trimmed.fq.gz"
    trim_report = "${sequence_id}_trim_report.txt"

    """
    cutadapt -a ${params.adapter_seq} -o ${trimmed_reads} -j ${params.trim_ncores} -q ${params.trim_qual} ${fastq_file} > ${trim_report}
    """

  stub:
    trimmed_reads = "${sequence_id}_trimmed.fq.gz"
    trim_report = "${sequence_id}_trim_report.json"

    """
    touch "${trimmed_reads}"
    touch "${trim_report}"
    """
}
    
