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

process basic_bwa {
  conda params.HOME_REPO + '/nf/envs/bwa.yaml'

  input
    tuple val(sequence_id), file(fq1), file(fq2), file(fq3)

  output:
    tuple val(sequence_id), file(bam)

  script:
    bam="${sequence_id}.bam"
    int_bam="${sequence_id}.uns.bam"
    """
    if [ "${fq2}" == "null" ]; then
        bwa mem -R "@RG\tID:${sequence_id}.RG\tSM:${sequence_id}\tPL:UNK\tPU:UNK}" \
                -t $params.alignment_ncore -k 17 $params.genome_reference $fq1 | samtools view -hb > ${int_bam}
    elif [ "${fq3}" == "null" ]; then
        bwa mem -R "@RG\tID:${sequence_id}.RG\tSM:${sequence_id}\tPL:UNK\tPU:UNK" \
                -t $params.alignment_ncore -k 17 $params.genome_reference $fq1 $fq2 | samtools view -hb > ${int_bam}
    else
        bwa mem -R "@RG\tID:${sequence_id}.RG\tSM:${sequence_id}\tPL:UNK\tPU:UNK" \
                -t $params.alignment_ncore -k 17 $params.genome_reference $fq1 $fq2 | samtools view -hb > ${sequence_id}.p.bam

        bwa mem -R "@RG\tID:${sequence_id}.RG\tSM:${sequence_id}\tPL:UNK\tPU:UNK" \
                -t $params.alignment_ncore -k 17 $params.genome_reference $fq3 | samtools view -hb > ${sequence_id}.w.bam

        samtools merge -o $int_bam ${sequence_id}.p.bam ${sequence_id}.w.bam
    fi
        

    samtools sort -@ $params.alignment_ncore "${int_bam}" -o "${bam}"
    """

   stub:
    bam="${sequence_id}.bam"
    int_bam="${sequence_id}.uns.bam"
    """
    touch $bam
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
    samtools stats "${bam_file}" | egrep "reads mapped:|reads unmapped:|reads MQ0:" | cut -f2,3 > stats_tmp1
    samtools view -h -q 20 "${bam_file}" | samtools stats | grep "reads mapped:" | sed 's/mapped:/mapped_Q20:/' | cut -f2,3 | awk '{print \$1"_"\$2,\$3}' >> stats_tmp1
    ok_cell=\$(samtools view -q 20 "${bam_file}" | cut -f1 | tr '|' '\t' | cut -f2 | tr ':' '\t' | awk '{print \$2,\$3}' | sort | uniq -c | awk '\$1 >= 5000' | wc -l)
    echo "N cells >= 5000 reads: \${ok_cell}" >> stats_tmp1
    good_cell=\$(samtools view -q 20 "${bam_file}" | cut -f1 | tr '|' '\t' | cut -f2 | tr ':' '\t' | awk '{print \$2,\$3}' | sort | uniq -c | awk '\$1 >= 20000' | wc -l)
    echo "N cells >= 20000 reads: \${good_cell}" >> stats_tmp1
    while read statline; do
      statname=\$(echo \$statline | tr ':' '\t' | cut -f1 | sed 's/^\\s\\+//g')
      statval=\$(echo \$statline | tr ':' '\t' | cut -f2 | sed 's/^\\s\\+//g')
      echo "${sequence_id},${bamname},${seqtype},${assay},${antibody},\${statname},\${statval}" >> "${alignment_stats}"
    done < stats_tmp1
    rm stats_tmp1
    """

  stub:
    bamname = bam_file.simpleName
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
    merged_file = base_name + ".alignment_qc.txt"
    """
    cmd=\$(echo "cat ${qc_files} @ ${merged_file}" | tr '@' '>')
    eval \$cmd
    """

  stub:
    merged_file = base_name + ".alignment_qc.txt"
    """
    touch "${merged_file}"
    """
} 

/*
 * This process converts a .bam file into a Signac-compatible
 * fragments file, ensures that it is sorted, block-compresses
 * and indexes the fragments file.
 *
 * Config-defined parameters
 * --------------------------
 * HOME_REPO - the path to the home repository
 * fragment_ncore - the number of cores to use
 */
process bam_to_frag {
  conda params.HOME_REPO + '/nf/envs/bamtofrag.yaml'
  input:
    file bam_file

  output:
    file fragment_file
    file fragment_index

  script:
    fragment_file = (bam_file.toString() - '.bam') + '.frag.tsv.gz'
    fragment_index = fragment_file + '.tbi'
    """
    samtools index $bam_file
    python $params.HOME_REPO/py/bam2frag.py $bam_file $fragment_file --ncores $params.fragment_ncore
    zcat $fragment_file | bedtools sort -i /dev/stdin | bgzip -c > ${fragment_file}.tmp
    mv ${fragment_file}.tmp ${fragment_file}
    tabix -p bed $fragment_file
    """

   stub:
     fragment_file = (bam_file.toString() - '.bam') + '.frag.tsv.gz'
     fragment_index = fragment_file + '.tbi'
     """
     touch $fragment_file $fragment_index
     """
}


/*
 * This process merges multiple fragment files
 *
 * Config-defined parameters
 * ------------
 * HOME_REPO - the path to the home repository
 */
process merge_frag_files {
  conda $params.HOME_REPO + '/nf/envs/bamtofrag.yaml'

  input:
    file fragment_file  // .collect() has been run
    val run_name

  output:
    file merged_fragments
    file merged_fragments_index

  script:
    merged_fragments = run_name + '_fragments_allAntibodies.tsv.gz'
    merged_fragments_idnex = merged_fragments + '.tbi'
    """
    zcat $fragment_file | bedtools sort - /dev/stin | bgzip -c > ${merged_fragments}
    tabix -p bed ${merged_fragments}
    """

  stub:
    merged_fragments = run_name + '_fragments_allAntibodies.tsv.gz'
    merged_fragments_idnex = merged_fragments + '.tbi'
    """
    touch $merged_fragments
    touch $merged_fragments_index
    """
}
  
   
