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
      tuple val(sequence_id), file(fastq_trimmed1), val(seqtype)
      

    output:
      tuple val(sequence_id), file(outbam), val(seqtype), val(assay), val(antibody)

  script: 
    prefix = fastq_trimmed1.simpleName
    outbam = prefix + 'Aligned.sortedByCoord.out.bam'
    input_fq1 = "${fastq_trimmed1}"
    assay = input_fq1.split("__")[1]
    antibody = input_fq1.split("__")[2]

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
    input_fq1 = "${fastq_trimmed1}"
    assay = input_fq1.split("__")[1]
    antibody = input_fq1.split("__")[2]
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
    tuple val(sequence_id), file(fq_file), val(seqtype)

  output:
    tuple val(sequence_id), file(aln_bam), val(seqtype), val(assay), val(antibody)

  script:
    fq_pfx = fq_file.simpleName
    // alignment_id = "${sequence_id}_${fq_pfx}"
    aln_bam = "${fq_pfx}.bam"
    seqid = fq_pfx.split("__")[0]
    assay = fq_pfx.split("__")[1]
    antibody = fq_pfx.split("__")[2]
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
    seqid = fq_pfx.split("__")[0]
    assay = fq_pfx.split("__")[1]
    antibody = fq_pfx.split("__")[2]
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


/*
 * This process extracts basic alignment QC metrics
 */
process alignment_qc {
  conda params.HOME_REPO + '/nf/envs/bwa.yaml'

  input:
    tuple val(sequence_id), file(bam_file), val(seqtype), val(assay), val(antibody)

  output:
    tuple val(sequence_id), file(alignment_stats)

  script:
    bamname = bam_file.simpleName
    alignment_stats = bamname - '.bam' + '.alignment_stats.txt'
    """
    samtools stats "${bam_file}" | egrep "reads mapped:|reads unmapped:|reads MQ0:" | cut -f2 > stats_tmp1
    samtools view -h -q 20 "${bam_file}" | samtools stats | grep "reads mapped:" | sed 's/mapped:/mapped Q20:/' | cut -f2 >> stats_tmp1
    ok_cell=$(samtools view -q 20 "${bam_file}" | cut -f1 | tr '|' '\t' | cut -f2 | tr ':' '\t' | awk '{print $2,$3}' | sort | uniq -c | awk '$1 >= 5000' | wc -l)
    echo "N cells >= 5000 reads: ${ok_cell}" >> stats_tmp1
    good_cell=$(samtools view -q 20 "${bam_file}" | cut -f1 | tr '|' '\t' | cut -f2 | tr ':' '\t' | awk '{print $2,$3}' | sort | uniq -c | awk '$1 >= 20000' | wc -l)
    echo "N cells >= 20000 reads: ${good_cell}" >> stats_tmp1
    while read statline; do
      statname=$(echo $statline | tr ':' '\t' | cut -f1 | sed 's/^\s\+//g')
      statval=$(echo $statline | tr ':' '\t' | cut -f2 | sed 's/^\s\+//g')
      echo "${sequence_id},${bamname},${seqtype},${assay},${antibody},${statname},${statval}" >> "${alignment_stats}"
    done < stats_tmp1
    rm stats_tmp1
    """

  stub:
    bamname = bam_file.simpelName
    alignment_stats = bamname - '.bam' + '.alignment_stats.txt'
    """
    touch "${alignment_stats}"
    """
}

/*
 * This process concatenates multiple alignment QC files into a single output file
 *
 * to be used as merge_alignment_qc(qc_files.collect())
 */
process merge_alignment_qc {
  input:
    file qc_files 
    val base_name

  output:
    file merged_file

  script:
    merged_file = base_name + '.alignment_qc.txt"
    """
    cat "${qc_files}" > "${merged_file}"
    """

  stub:
    merged_file = base_name + '.alignment_qc.txt"
    """
    touch "${merged_file}"
    """
} 
