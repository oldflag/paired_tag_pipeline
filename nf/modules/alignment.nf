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
    sample = input_fq1.split("__")[3]
    tmp_fq = "trimmed.fq"
    trim_report = prefix + ".trim_report.txt"
    """
    mkfifo ${tmp_fq}
    cutadapt -a ${params.adapter_seq} -g ${params.universal_seq} -g ${params.transposase_seq} --times 2 -o ${tmp_fq} -j ${params.trim_ncores} -q ${params.trim_qual} -m 25 ${input_fq1} > ${trim_report} &
    STAR --readFilesIn $tmp_fq \\
        --runThreadN $params.alignment_ncore \\
        --twopassMode None \\
        --outSAMtype BAM SortedByCoordinate \\
        --readFilesCommand cat \\
        --outSAMunmapped Within \\
        --limitBAMsortRAM $params.ramsize \\
        --genomeDir $params.star_index \\
        --outFileNamePrefix "${prefix}"
    rm ${tmp_fq}
    """

  stub:
    prefix = fastq_trimmed1.simpleName
    outbam = prefix + 'Aligned.sortedByCoord.out.bam'
    input_fq1 = "${fastq_trimmed1}"
    assay = input_fq1.split("__")[1]
    antibody = input_fq1.split("__")[2]
    sample = input_fq1.split("__")[3]
    trim_report = prefix + ".trim_report.txt"
    """
    touch "${outbam}" "${trim_report}"
    """
}


process basic_bwa {
  conda params.HOME_REPO + '/nf/envs/bwa.yaml'

  input:
    tuple val(sequence_id), val(fq1), val(fq2), val(fq3)

  output:
    tuple val(sequence_id), file(bam)

  script:
    bam="${sequence_id}.bam"
    int_bam="${sequence_id}.uns.bam"
    """
    if [ "${fq2}" == "null" ]; then
        bwa mem -R "@RG\\tID:${sequence_id}.RG\\tSM:${sequence_id}\\tPL:UNK\\tPU:UNK}" \
                -t $params.alignment_ncore -k 17 $params.genome_reference $fq1 | samtools view -hb > ${int_bam}
    elif [ "${fq3}" == "null" ]; then
        bwa mem -R "@RG\\tID:${sequence_id}.RG\\tSM:${sequence_id}\\tPL:UNK\\tPU:UNK" \
                -t $params.alignment_ncore -k 17 $params.genome_reference $fq1 $fq2 | samtools view -hb > ${int_bam}
    else
        bwa mem -R "@RG\\tID:${sequence_id}.RG\\tSM:${sequence_id}\\tPL:UNK\\tPU:UNK" \
                -t $params.alignment_ncore -k 17 $params.genome_reference $fq1 $fq2 | samtools view -hb > ${sequence_id}.p.bam

        bwa mem -R "@RG\\tID:${sequence_id}.RG\\tSM:${sequence_id}\\tPL:UNK\\tPU:UNK" \
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


process trimming_bwa {
  conda params.HOME_REPO + '/nf/envs/bwa.yaml'

  input:
    tuple val(sequence_id), val(fq1), val(fq2), val(fq3)

  output:
    tuple val(sequence_id), file(bam)

  script:
    bam="${sequence_id}.bam"
    int_bam="${sequence_id}.uns.bam"
    trim_report = "${sequence_id}.trim_report.txt"
    """
    if [ "${fq2}" == "null" ]; then
        tmp_fq="trimmed.fq"
        mkfifo \$tmp_fq
        cutadapt -a ${params.adapter_seq} -g ${params.universal_seq} -g ${params.transposase_seq} --times 2 -o ${fq1} -j ${params.trim_ncores} -q ${params.trim_qual} -m 25 ${fq1} > ${trim_report} &
        bwa mem -R "@RG\\tID:${sequence_id}.RG\\tSM:${sequence_id}\\tPL:UNK\\tPU:UNK}" \\
                -t $params.alignment_ncore -k 17 $params.genome_reference \$tmp_fq | samtools sort -@ $params.alignment_ncore -o "${bam}"
        rm \$tmp_fq
    elif [ "${fq3}" == "null" ]; then
        tmp_fq1="tmp_r1.fq"
        tmp_fq2="tmp_r2.fq"
        cutadapt -a ${params.adapter_seq} -g ${params.universal_seq} -g ${params.transposase_seq} --times 2 -o \${tmp_fq1} -p \${tmp_fq2} -j ${params.trim_ncores} -q ${params.trim_qual} -m 25 ${fq1} ${fq2} > ${trim_report} 
        bwa mem -R "@RG\\tID:${sequence_id}.RG\\tSM:${sequence_id}\\tPL:UNK\\tPU:UNK" \\
                -t $params.alignment_ncore -k 17 $params.genome_reference \$tmp_fq1 \$tmp_fq2 | samtools sort -@ $params.alignment_ncore -o "${bam}"
        rm \$tmp_fq1 \$tmp_fq2
    else
        tmp_fq1="tmp_r1.fq"
        tmp_fq2="tmp_r2.fq"
        tmp_wid="tmp_wid.fq"
        cutadapt -a ${params.adapter_seq} -o \${tmp_wid} -j ${params.trim_ncores} -q ${params.trim_qual} -m 25 ${fq3} > trim_widow.text 
        cutadapt -a ${params.adapter_seq} -j ${params.trim_ncores} -q ${params.trim_qual} -o \$tmp_fq1 -p \$tmp_fq2 $fq1 $fq2 > trim_report.txt 

        bwa mem -R "@RG\\tID:${sequence_id}.RG\\tSM:${sequence_id}\\tPL:UNK\\tPU:UNK" \
                -t $params.alignment_ncore -k 17 $params.genome_reference \$tmp_fq1 \$tmp_fq2 | samtools view -hb > ${sequence_id}.p.bam

        bwa mem -R "@RG\\tID:${sequence_id}.RG\\tSM:${sequence_id}\\tPL:UNK\\tPU:UNK" \
                -t $params.alignment_ncore -k 17 $params.genome_reference \$tmp_wid | samtools view -hb > ${sequence_id}.w.bam

        samtools merge -o $int_bam ${sequence_id}.p.bam ${sequence_id}.w.bam
        rm \$tmp_fq1 \$tmp_fq2 \$tmp_wid ${sequence_id}.p.bam \${sequence_id}.w.bam
    fi

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
    assay = fq_pfx.split("__")[1]
    antibody = fq_pfx.split("__")[2]
    sample = fq_pfx.split("__")[3]
    // alignment_id = "${sequence_id}_${fq_pfx}"
    aln_bam = "${fq_pfx}.bam"
    trim_report = "${fq_pfx}.trim_report.txt"
    tmp_fq = "temp.fq"
    """
    mkfifo ${tmp_fq}
    cutadapt -a ${params.adapter_seq} -o ${tmp_fq} -j ${params.trim_ncores} -q ${params.trim_qual} -m 25 ${fq_file} > ${trim_report} &
    bwa mem -t "${params.alignment_ncore}" "${params.genome_reference}" "${tmp_fq}" | samtools view -Shu - | samtools sort - > "${aln_bam}"
    rm ${tmp_fq}
    d=\$(samtools view "${aln_bam}" | head -n 10 | wc -l)
    if [[ "\${d}" -lt 2 ]]; then
      exit 255
    fi
    """

  stub:
    fq_pfx = fq_file.simpleName
    assay = fq_pfx.split("__")[1]
    antibody = fq_pfx.split("__")[2]
    sample = fq_pfx.split("__")[3]
    // alignment_id = "${sequence_id}_${fq_pfx}"
    aln_bam = "${fq_pfx}.bam"
    trim_report = "${fq_pfx}.trim_report.txt"
    """
    touch "${aln_bam}" "${trim_report}"
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
    n_bc=\$(samtools view "${bam_file}" | cut -f1 | tr '|' '\t' | cut -f2 | tr ':' '\t' | cut -f2,3 | sort | uniq | wc -l)
    aln_bc=\$(samtools view -F 4 "${bam_file}" | cut -f1 | tr '|' '\t' | cut -f2 | tr ':' '\t' | cut -f2,3 | sort | uniq | wc -l)
    aln_q20_bc=\$(samtools view -q 20 "${bam_file}" | cut -f1 | tr '|' '\t' | cut -f2 | tr ':' '\t' | cut -f2,3 | sort | uniq | wc -l)
    echo "N barcodes: \${n_bc}" >> stats_tmp1
    echo "N aligned barcodes: \${aln_bc}" >> stats_tmp1
    echo "N barcodes Q20: \${aln_q20_bc}" >> stats_tmp1
    ok_cell=\$(samtools view -q 20 "${bam_file}" | cut -f1 | tr '|' '\t' | cut -f2 | tr ':' '\t' | awk '{print \$2,\$3}' | sort | uniq -c | awk '\$1 >= 5000' | wc -l)
    mean_aln_cell=\$(samtools view -F 4 "${bam_file}" | cut -f1 | tr '|' '\t' | cut -f2 | tr ':' '\t' | awk '{print \$2, \$3}' | sort | uniq -c | awk 'BEGIN{s=0; n=1}{s = s + \$1; n=n+1}END{print(s/n)}')
    mean_cell_q20=\$(samtools view -q 20 "${bam_file}" | cut -f1 | tr '|' '\t' | cut -f2 | tr ':' '\t' | awk '{print \$2, \$3}' | sort | uniq -c | awk 'BEGIN{s=0; n=1}{s = s + \$1; n=n+1}END{print(s/n)}')
    echo "Mean aligned per barcode: \${mean_aln_cell}" >> stats_tmp1
    echo "Mean Q20 aligned per barcode: \${mean_cell_q20}" >> stats_tmp1
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
    file fragment_log
    file wig_track
    file fragsize_hist

  script:
    fragment_file = (bam_file.toString() - '.bam') + '.frag.tsv.gz'
    fragment_index = fragment_file + '.tbi'
    fragment_log = fragment_file + '.log'
    fragsize_hist = (bam_file.toString() - '.bam') + '.fragsize_hist.txt'
    wig_track = (bam_file.toString() - '.bam') + '.bw'
    """
    samtools index $bam_file
    samtools view -H $bam_file | grep @SQ | sed 's/@SQ\tSN://g' | sed 's/LN://g' > genome.txt
    python $params.HOME_REPO/py/bam2frag.py $bam_file $fragment_file --ncores $params.fragment_ncore --nocb | tee $fragment_log
    zcat $fragment_file | bedtools sort -i /dev/stdin | bgzip -c > ${fragment_file}.tmp
    mv ${fragment_file}.tmp ${fragment_file}
    zcat "${fragment_file}" | awk '{print \$1"\t"\$2"\t"\$3"\t"\$4"\t.\t."}' > unsorted.bed
    bedtools sort -i unsorted.bed > sorted.bed
    bedtools genomecov -i sorted.bed -g genome.txt -bg > coverage.bg
    bedGraphToBigWig coverage.bg genome.txt $wig_track
    samtools view -f 3 "${bam_file}" | cut -f1,9 | awk '\$2 > 0 && \$2 < 2000' | tr '|' '\t' | awk '{print \$(NF-1)"\t"\$(NF)}' | sort | bash "${params.HOME_REPO}/sh/average.sh" /dev/stdin | cut -f2 | sort | uniq -c | awk '{print \$2"\t"\$1}' | sort -n -k1,1 > frags.tmp
    cat frags.tmp | awk -v bam="${bam_file}" '{print bam"\t"\$1"\t"\$2}' > "${fragsize_hist}"
    """

   stub:
    fragment_file = (bam_file.toString() - '.bam') + '.frag.tsv.gz'
    fragment_index = fragment_file + '.tbi'
    fragment_log = fragment_file + '.log'
    wig_track = (bam_file.toString() - '.bam') + '.bw'
    fragsize_hist = (bam_file.toString() - '.bam') + '.fragsize_hist.txt'
    """
    touch $fragment_file $fragment_log $fragment_index $wig_track $fragsize_hist
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
    merged_fragments_index = merged_fragments + '.tbi'
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
  
   
