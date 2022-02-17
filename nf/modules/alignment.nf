/*
 * Modules of STAR(RNA) alignment for single and paired reads
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
 *   + ramsize: size of RAM
 */
process star_aligner_single {
  conda params.HOME_REPO + '/nf/envs/star.yaml'
  
    input:
      tuple val(sequence_id), file(fastq_trimmed1), val(keyvalues)
      

    output:
      tuple val(sequence_id), file(outbam), val(keyvalues)

  script: 
    prefix = fastq_trimmed1.simpleName
    outbam = prefix + 'Aligned.sortedByCoord.out.bam'
    input_fq1 = "${fastq_trimmed1}"

    """
    STAR --readFilesIn $input_fq1 \\
        --runThreadN $params.alignment_ncore \\
        --twopassMode Basic \\
        --outSAMtype BAM SortedByCoordinate \\
        --readFilesCommand zcat \\
        --outSAMunmapped Within \\
        --limitBAMsortRAM $params.ramsize \\
        --genomeDir $params.star_index \\
        --outFileNamePrefix "${prefix}"

    """

  stub:
    prefix = fastq_trimmed1.simpleName
    outbam = prefix + 'Aligned.sortedByCoord.out.bam'
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
    tuple val(sequence_id), file(fq_file), val(keyvalues)

  output:
    tuple val(sequence_id), file(aln_bam), val(keyvalues)

  script:
    fq_pfx = fq_file.simpleName
    // alignment_id = "${sequence_id}_${fq_pfx}"
    aln_bam = "${fq_pfx}.bam"
    """
    bwa mem -t "${params.alignment_ncore}" "${params.genome_reference}" "${fq_file}" | samtools view -hb > "${aln_bam}"
    d=\$(samtools view "${aln_bam}" | head -n 10 | wc -l)
    if [[ "\${d}" -lt 2 ]]; then
      exit 255
    fi
    """

  stub:
    fq_pfx = fq_file.simpleName
    // alignment_id = "${sequence_id}_${fq_pfx}"
    aln_bam = "${fq_pfx}.bam"
    """
    touch "${aln_bam}"
    """
}


/*
 * This process defines how to merge multiple bam files together. It
 * is written to be used with a .collect() statement; i.e.,
 *
 * mbam = merge_bams(sorted_bams.collect())
 *
 */
process merge_bams {
  
  conda params.HOME_REPO + '/nf/envs/bwa.yaml'

  input:
    tuple file(bam_file), val(base_name1), val(base_name2)

  output:
    tuple file(merged_bam), val(base_name1), val(base_name2)

  script:
    bamlist = bam_file.join('\n')
    base_name = base_name1+base_name2
    merged_bam = base_name+ '.bam'
    """
    echo "${bamlist}" > bamlist.txt
    samtools merge -b bamlist.txt "${merged_bam}"
    """

  stub:
    bamlist = bam_file.join('\n')
    base_name = base_name1+base_name2
    merged_bam = base_name + '.bam'
    """
    touch "${merged_bam}"
    """
}
