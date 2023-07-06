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
  // conda params.HOME_REPO + '/nf/envs/star.yaml'
  
    input:
      tuple val(sequence_id), file(fastq_trimmed1), val(seqtype)
      file star_index
      

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
        --genomeDir "${star_index}" \\
        --outFileNamePrefix "${prefix}" \\
        --genomeLoad LoadAndKeep
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
  // conda params.HOME_REPO + '/nf/envs/bwa.yaml'

  input:
    tuple val(sequence_id), file(fq_file), val(seqtype)
    file genome_reference
    file bwa_index

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
    //"${bwa_index}/${genome_reference}"
    
    """
    mkfifo ${tmp_fq}
    cutadapt -a ${params.adapter_seq} -o ${tmp_fq} -j ${params.trim_ncores} -q ${params.trim_qual} -m 25 ${fq_file} > ${trim_report} &
    bwa mem -t "${params.alignment_ncore}" "${bwa_index}/${genome_reference}" "${tmp_fq}" | samtools view -Shu - | samtools sort - > "${aln_bam}"
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
  
  // conda params.HOME_REPO + '/nf/envs/bwa.yaml'

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
  // conda params.HOME_REPO + '/nf/envs/bwa.yaml'

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

process add_tags {
  // conda params.HOME_REPO + '/nf/envs/bwa.yaml'
  
  input:
    tuple val(alignment_id), file(bam_file), val(seqtype), val(assay_id), val(antibody)
    tuple file(annot1_file), val(annot1_tag), val(annot1_fmt)
    tuple file(annot2_file), val(annot2_tag), val(annot2_fmt)
    tuple file(annot3_file), val(annot3_tag), val(annot3_fmt)
    tuple file(annot4_file), val(annot4_tag), val(annot4_fmt)
    tuple file(annot5_file), val(annot5_tag), val(annot5_fmt)
    tuple file(annot6_file), val(annot6_tag), val(annot6_fmt)
    file py_dir


  output:
    tuple val(alignment_id), file(annot_bam_file), val(seqtype), val(assay_id), val(antibody)

  script:
    annot_bam_file = bam_file.simpleName - '.bam' + '_tag.bam'
    sample_id=bam_file.simpleName.split('__')[3]
    """
    strip () {
        echo "\${1%%.*}"
    }

    tag_args=""
    if [  \$(strip "${annot1_file}") != "input" ]; then
        tag_args="\${tag_args} --track ${annot1_file}:${annot1_tag}:${annot1_fmt}"
    fi
    
    if [  \$(strip "${annot2_file}") != "input" ]; then
        tag_args="\${tag_args} --track ${annot2_file}:${annot2_tag}:${annot2_fmt}"
    fi
    if [  \$(strip "${annot3_file}") != "input" ]; then
        tag_args="\${tag_args} --track ${annot3_file}:${annot3_tag}:${annot3_fmt}"
    fi
    if [  \$(strip "${annot4_file}") != "input" ]; then
        tag_args="\${tag_args} --track ${annot4_file}:${annot4_tag}:${annot4_fmt}"
    fi
    if [  \$(strip "${annot5_file}") != "input" ]; then
        tag_args="\${tag_args} --track ${annot5_file}:${annot5_tag}:${annot5_fmt}"
    fi
    if [  \$(strip "${annot6_file}") != "input" ]; then
        tag_args="\${tag_args} --track ${annot6_file}:${annot6_tag}:${annot6_fmt}"
    fi

    if [ "\${tag_args}" == "" ]; then
        python "${py_dir}"/add_tags.py "${bam_file}" "${annot_bam_file}" --library "${alignment_id}" --antibody "${antibody}" --assay_id "${assay_id}" --sample_id "${sample_id}"
    else
        python "${py_dir}"/add_tags.py "${bam_file}" "${annot_bam_file}" --library "${alignment_id}" --antibody "${antibody}" --assay_id "${assay_id}" --sample_id "${sample_id}" \$tag_args
    fi
    """

  stub:
    annot_bam_file = bam_file.simpleName.replace('.bam', '_tag.bam')
    """
    touch "${annot_bam_file}"
    """
}

