/*
 * Modules of Star alignment for single and paired reads
 */

nextflow.enable.dsl=2

/*
 * This module defines a process for STAR alignment
 *
 * Config-defined parameters:
 * --------------------------
 *   + fastq_trimmed1: read sequence
 *   + fastq_trimmed2: read sequence
 *   + star_index: index
 *   + alignment_ncore: the number of cores to use for alignment
 */
process star_aligner_single {
  conda params.HOME_REPO + '/nf/envs/star.yaml'
  
    input:
      tuple val(sequence_id), val(fastq_trimmed1)
      

    output:
      tuple val(sequence_id), path ('*.bam')

  script: 
    prefix = "${sequence_id}"
    prefix = "${sequence_id}"
    outbam = prefix + '.Aligned.sortedByCoord.out.bam'
    input_fq1 = params.datadir + "${fastq_trimmed1}"

    """
    STAR --readFilesIn $input_fq1 \\
        --runThreadN $params.alignment_ncore \\
        --twopassMode Basic \\
        --outSAMtype BAM SortedByCoordinate \\
        --readFilesCommand zcat \\
        --outSAMunmapped Within \\
        --limitBAMsortRAM 16000000000 \\
        --genomeDir $params.star_index \\
        --outFileNamePrefix $bamPrefix

    """

  stub:
    prefix = "${sequence_id}"
    outbam = prefix + '.Aligned.sortedByCoord.out.bam'
    """
    touch "${outbam}"
    """
}

process star_aligner_pair {
  conda params.HOME_REPO + '/envs/star.yaml'
  
    input:
      tuple val(sequence_id), val(fastq_trimmed1), val(fastq_trimmed2)
      

    output:
      path ('*.bam'), emit: star_aligned
      path "*.out", emit: alignment_report
      file "*SJ.out.tab"
      file "*Log.out"
      file "*.bam.bai"
      file "*.txt"


  script: 
    prefix = "${sequence_id}"
    bamPrefix =  params.outputdir + prefix
    input_fq1 = params.datadir + "${fastq_trimmed1}"
    input_fq2 = params.datadir + "${fastq_trimmed2}"

    """
    STAR --readFilesIn $input_fq1 $input_fq2 \\
        --runThreadN "${params.alignment_ncore}" \\
        --twopassMode Basic \\
        --outSAMtype BAM SortedByCoordinate \\
        --readFilesCommand zcat \\
        --outSAMunmapped Within \\
        --limitBAMsortRAM 16000000000 \\
        --genomeDir $params.star_index \\
        --outFileNamePrefix $bamPrefix

    """

  stub:
    prefix = "${sequence_id}"
    bamPrefix =  params.outputdir + prefix
    outbam = bamPrefix + '.Aligned.sortedByCoord.out.bam'
    """
    touch "${outbam}"
    """
}

/*
 * This process defines a mechanism for aligning genomic dna
 * with bwa
 *
 * Config-defined parameters:
 * ---------------------------
 * alignment_ncore - the number of cores to use in alignment
 * genome_reference - the reference file to use for bwa alignment
 *
 * TODO: parameter search
 */

process bwa_aligner_single {
  conda params.HOME_REPO + '/nf/envs/bwa.yaml'

  input:
    tuple val(sequence_id), file(fq_file)

  output:
    tuple val(sequence_id), file(aln_bam)

  script:
    aln_bam = "${sequence_id}.bam"
    rbase = params.genome_reference.toString().strip('.gz').strip('.fa')
    """
    bwa mem -t "${params.alignment_ncore}" "${rbase}" "${fq_file}" | samtools view -hb > "${aln_bam}"
    """

  stub:
    aln_bam = "${sequence_id}".bam
    """
    touch "${aln_bam}"
    """
}
