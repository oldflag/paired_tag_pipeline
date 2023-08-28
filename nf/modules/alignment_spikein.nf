/*
 * Modules of STAR(RNA) alignment for single and paired reads
 */

nextflow.enable.dsl=2

def checkExists(var_) {
  try {
    varName.toString()
    return true
  } catch(ex) {
    return false
  }
}

/*
 * This module defines a process for STAR alignment given
 * both the primary and spike-in species (hs/mm)
 *
 */
process star_aligner {
  // conda params.HOME_REPO + '/nf/envs/star.yaml'
  
    input:
      tuple val(sequence_id), file(fastq_trimmed1), file(fastq_trimmed2), val(seqtype)
      path index_primary, stageAs:'temp/*'
      path index_spike, stageAs:'temp2/*'
      file py_dir
      

    output:
      tuple val(sequence_id), file(outbam), val(seqtype), val(assay), val(antibody)
      tuple val(sequence_id), file(outbam_spike), val(seqtype), val(assay), val(antibody)
      tuple val(sequence_id), file(spikein_info)
      
   

  script: 
    prefix = fastq_trimmed1.simpleName
    // star suffix: Aligned.sortedByCoord.out.bam
    outbam = prefix + 'Aligned.sortedByCoord.out.bam'
    outbam_spike = prefix + '_spikeAligned.sortedByCoord.out.bam'
    input_fq1 = "${fastq_trimmed1}"
    assay = input_fq1.split("__")[1]
    antibody = input_fq1.split("__")[2]
    sample = input_fq1.split("__")[3]
    tmp_fq = "trimmed.fq"
    // index_primary = params.star_index[params.SPECIES]
    // index_spike = params.star_index[params.SPIKEIN_SPECIES]
    spikein_info = "${sequence_id}.${assay}.${antibody}.${sample}.${seqtype}.spikein_info.txt"
    trim_report = prefix + ".trim_report.txt"
    if ( checkExists(params.NO_SPIKEIN) && params.NO_SPIKEIN == "yes" ) {
        ext_args=" --duplicate_reads"
    } else {
        ext_args=""
    }
    """
    cutadapt -a ${params.adapter_seq} -g ${params.universal_seq} -g ${params.transposase_seq} --times 2 -o ${tmp_fq} -j ${params.trim_ncores} -q ${params.trim_qual} -m 20 ${input_fq1} > ${trim_report}
    #python "${py_dir}/filt_missing.py" $tmp_fq $fastq_trimmed2 r1_in.fq r2_in.fq
    cat ${tmp_fq} | sed 's/^\$/N/g' | awk 'BEGIN { RS = "@"; ORS = "" } index(\$1, ":rna|") { print "@" \$0 }' | gzip -c > foo && mv foo "${tmp_fq}.gz"
    #rm $tmp_fq
    STAR --readFilesIn "${tmp_fq}.gz" \\
        --runThreadN $params.alignment_ncore \\
        --twopassMode None \\
        --outSAMtype BAM Unsorted \\
        --readFilesCommand zcat \\
        --outSAMunmapped Within KeepPairs \\
        --limitBAMsortRAM $params.ramsize \\
        --genomeDir $index_primary \\
        --outSAMmultNmax 1 \\
        --outFileNamePrefix "${prefix}_uns"
    STAR --readFilesIn "${tmp_fq}.gz" \\
        --runThreadN $params.alignment_ncore \\
        --twopassMode None \\
        --outSAMtype BAM Unsorted \\
        --readFilesCommand zcat \\
        --outSAMunmapped Within KeepPairs \\
        --limitBAMsortRAM $params.ramsize \\
        --genomeDir $index_spike \\
        --outSAMmultNmax 1 \\
        --outFileNamePrefix "${prefix}_spike_uns"
    # star suffix: Aligned.sortedByCoord.out.bam
    #rm r1_in.fq r2_in.fq
    samtools sort -n "${prefix}_unsAligned.out.bam" -o "${prefix}_nameSorted.bam"
    samtools sort -n "${prefix}_spike_unsAligned.out.bam" -o "${prefix}_spike_nameSorted.bam"
    python "${py_dir}/assign_spikeins.py" \\
         "${prefix}_nameSorted.bam" \\
         "${prefix}_spike_nameSorted.bam" \\
         "${prefix}_uns.bam" \\
         "${prefix}_spike_uns.bam" \\
         --info "${spikein_info}" \\
         --seqtype rna ${ext_args} \\
         --unpaired
    rm "${prefix}_unsAligned.out.bam" "${prefix}_spike_unsAligned.out.bam" "${prefix}_nameSorted.bam" "${prefix}_spike_nameSorted.bam"
    samtools sort "${prefix}_uns.bam" -o "${outbam}"
    samtools sort "${prefix}_spike_uns.bam" -o "${outbam_spike}"
    rm "${prefix}_uns.bam" "${prefix}_spike_uns.bam"
    
    """

  stub:
    prefix = fastq_trimmed1.simpleName
    outbam = prefix + 'Aligned.sortedByCoord.out.bam'
    outbam_spike = prefix + '_spikeAligned.sortedByCoord.out.bam'
    input_fq1 = "${fastq_trimmed1}"
    assay = input_fq1.split("__")[1]
    antibody = input_fq1.split("__")[2]
    sample = input_fq1.split("__")[3]
    trim_report = prefix + ".trim_report.txt"
    spikein_info = "${sequence_id}.${assay}.${antibody}.${sample}.${seqtype}.spikein_info.txt"
    """
    touch "${outbam}" "${trim_report}" "${outbam_spike}" "${spikein_info}" 
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

process bwa_aligner {
  // conda params.HOME_REPO + '/nf/envs/bwa.yaml'

  input:
    tuple val(sequence_id), file(fq_file), file(fq_file2),val(seqtype)
    path primary_bwa_index, stageAs:'temp/*'
    file primary_ref 
    path spike_bwa_index, stageAs:'temp2/*'
    file spike_ref
    file py_dir

  output:
    tuple val(sequence_id), file(aln_bam), val(seqtype), val(assay), val(antibody)
    tuple val(sequence_id), file(aln_spike_bam), val(seqtype), val(assay), val(antibody)
    tuple val(sequence_id), file(spikein_info)
    

  script:
    fq_pfx = fq_file.simpleName
    assay = fq_pfx.split("__")[1]
    antibody = fq_pfx.split("__")[2]
    sample = fq_pfx.split("__")[3]
    // alignment_id = "${sequence_id}_${fq_pfx}"
    aln_bam = "${fq_pfx}.bam"
    aln_spike_bam = "${fq_pfx}_spike.bam"
    init_bam = "${fq_pfx}.init.bam"
    init_spike_bam = "${fq_pfx}_spike.init.bam"
    filt_bam = "${fq_pfx}.filt.bam"
    filt_spike_bam = "${fq_pfx}_spike.filt.bam"
    trim_report = "${fq_pfx}.trim_report.txt"
    tmp_fq = "temp.fq"
    spikein_info = "${sequence_id}.${assay}.${seqtype}.spikein_info.txt"
    // primary_ref = params.genome_reference[params.SPECIES]
    // spike_ref = params.genome_reference[params.SPIKEIN_SPECIES]
    if ( checkExists(params.NO_SPIKEIN) && params.NO_SPIKEIN == "yes" ) {
        ext_args=" --duplicate_reads"
    } else {
        ext_args=""
    }
    """
    cutadapt -a ${params.adapter_seq} -o ${tmp_fq} -j ${params.trim_ncores} -q ${params.trim_qual} -m 0 ${fq_file} > ${trim_report} 
    cat ${tmp_fq} | sed 's/^\$/N/g' | awk 'BEGIN { RS = "@"; ORS = "" } index(\$1, ":dna|") { print "@" \$0 }' | gzip -c > foo && mv foo "${tmp_fq}.gz"
    alen=\$(zcat "${fq_file}" | awk 'NR % 4 == 2' | head -n 10 | awk 'BEGIN{tot=0}{tot += length(\$1)}END{print(tot/10)}')
    # if the average read length is long enough, try to align in paired-end mode
    if [ "\${alen}" -ge 150 ]; then
      zcat ${fq_file2} | sed 's/^\$/N/g' | awk 'BEGIN { RS = "@"; ORS = "" } index(\$1, ":dna|") { print "@" \$0 }' | gzip -c > foo && mv foo read2.fq.gz
      bwa mem -T 0 -h 1 -t "${params.alignment_ncore}" "${primary_bwa_index}/${primary_ref}" "${tmp_fq}.gz" "read2.fq.gz" | samtools view -Shu - > primary.bam
      bwa mem -T 0 -h 1 -t "${params.alignment_ncore}" "${spike_bwa_index}/${spike_ref}" "${tmp_fq}.gz" "read2.fq.gz" | samtools view -Shu - > spike.bam
      view_xargs="-f 64"
    else
      bwa mem -T 0 -h 1 -t "${params.alignment_ncore}" "${primary_bwa_index}/${primary_ref}" "${tmp_fq}.gz" | samtools view -Shu - > primary.bam
      bwa mem -T 0 -h 1 -t "${params.alignment_ncore}" "${spike_bwa_index}/${spike_ref}" "${tmp_fq}.gz" | samtools view -Shu - > spike.bam
      view_xargs=""
    fi
    rm "${tmp_fq}.gz"
    samtools view -F 2048 -hb primary.bam | samtools sort -n -o "${init_bam}"
    samtools view -F 2048 -hb spike.bam | samtools sort -n -o "${init_spike_bam}"
    python "${py_dir}/assign_spikeins.py" "${init_bam}" "${init_spike_bam}" "${filt_bam}" "${filt_spike_bam}" --info "${spikein_info}" --seqtype dna ${ext_args}
    rm ${init_bam} ${init_spike_bam} primary.bam spike.bam
    samtools sort "${filt_bam}" -o "${aln_bam}"
    samtools sort "${filt_spike_bam}" -o "${aln_spike_bam}"
    rm ${filt_bam} ${filt_spike_bam}
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
    aln_spike_bam = "${fq_pfx}_spike.bam"
    trim_report = "${fq_pfx}.trim_report.txt"
    spikein_info = "${sequence_id}.${assay}.${antibody}.${sample}.${seqtype}.spikein_info.txt"
    """
    touch "${aln_bam}" "${aln_spike_bam}" "${trim_report}" "${spikein_info}"
    """
}

/*
 * This process extracts basic alignment QC metrics
 */
process alignment_qc {
  // conda params.HOME_REPO + '/nf/envs/pysam.yaml'

  input:
    tuple val(sequence_id), file(bam_file), file(spikein_bam) val(seqtype), val(assay), val(antibody)

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
  // conda params.HOME_REPO + '/nf/envs/bamtofrag.yaml'
  input:
    file bam_file
    file py_dir

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
    fragsize_hist = (bam_file.toString() - '.bam') + '.fragsize_histogram.txt'
    pyfh = fragsize_hist - '.txt' + '.fragscript.txt'
    wig_track = (bam_file.toString() - '.bam') + '.bw'
    """
    samtools index $bam_file
    samtools view -H $bam_file | grep @SQ | sed 's/@SQ\tSN://g' | sed 's/LN://g' > genome.txt
    head genome.txt
    python ${py_dir}/bam2frag.py $bam_file $fragment_file --ncores $params.fragment_ncore --nocb --fragsize_hist "${pyfh}" | tee $fragment_log
    zcat $fragment_file | bedtools sort -i /dev/stdin | bgzip -c > ${fragment_file}.tmp
    mv ${fragment_file}.tmp ${fragment_file}
    zcat "${fragment_file}" | awk '{print \$1"\t"\$2"\t"\$3"\t"\$4"\t.\t."}' > unsorted.bed
    bedtools sort -i unsorted.bed > sorted.bed
    bedtools genomecov -i sorted.bed -g genome.txt -bg > coverage.bg
    bedGraphToBigWig coverage.bg genome.txt $wig_track
    samtools view -f 3 "${bam_file}" | cut -f1,9 | awk '\$2 > 0 && \$2 < 2000' | tr '|' '\t' | cut -f3,4 | sort | bash "${params.HOME_REPO}/sh/average.sh" /dev/stdin | cut -f2 | sort | uniq -c | awk '{print \$2"\t"\$1}' | sort -n -k1,1 > frags.tmp
    cat frags.tmp | awk -v bam="${bam_file}" '{print bam"\t"\$1"\t"\$2}' > "${fragsize_hist}"
    """

   stub:
     fragment_file = (bam_file.toString() - '.bam') + '.frag.tsv.gz'
     fragment_index = fragment_file + '.tbi'
     fragment_log = fragment_file + '.log'
     fragsize_hist = (bam_file.toString() - '.bam') + '.fragsize_histogram.txt'
     wig_track = (bam_file.toString() - '.bam') + '.bw'
     """
     touch $fragment_file $fragment_log $fragment_index $wig_track $fragsize_hist
     """
}


