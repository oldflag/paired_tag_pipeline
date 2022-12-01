/*
 * Modules for peak calling
 */
nextflow.enable.dsl=2

/*
 * This defines a process for peak calling from chip-seq or atac
 * bam files (dna)
 *
 * Config-defined parameters:
 * ----------------------------
 * HOME_REPO - location of the repository
 * macs_genome_type - the genome type (hs - human / mm - mouse)
 *
 */
process MACS2_peakcall {
  // conda params.HOME_REPO  + '/nf/envs/macs2.yaml'

  input:
    tuple file(bam_file), val(experiment_name), val(antibody_name)

  output:
    tuple val(experiment_name), file(merged_peaks), file(peaks_saf), val(antibody_name)

  script:
    basename = "${experiment_name}_${antibody_name}"
    narrow_peaks = "${basename}_peaks.narrowPeak"
    broad_peaks = "${basename}_peaks.broadPeak"
    merged_peaks = "${basename}_peaks_merged.bed"
    peaks_saf = "${basename}_peaks_merged.saf"
    filt_bam = "${basename}_filt.bam"
    """
    samtools view -hb -q 30 "${bam_file}" > "${filt_bam}" 
    macs2 callpeak -t "${filt_bam}" -n "${basename}" --outdir . \
       -q 0.1 -g "${params.macs_genome_type}" --nomodel
    macs2 callpeak -t "${filt_bam}" --broad -n "${basename}" --outdir . \
       -q 0.1 -g "${params.macs_genome_type}" --nomodel

    cat "${narrow_peaks}" "${broad_peaks}" | cut -f1-5 | bedtools sort -i /dev/stdin | bedtools merge -i /dev/stdin > "${merged_peaks}"
    echo "GeneID	Chr	Start	End	Strand" > "${peaks_saf}"
    cat "${merged_peaks}" | awk '{print "${antibody_name}_Peak."NR"@"\$1":"\$2":"\$3"\t"\$1"\t"\$2"\t"\$3"\t."}' >> "${peaks_saf}"
      nl=\$(cat ${peaks_saf} | wc -l)
      if [ "\${nl}" -eq "1" ]; then
          # just the header
          if [ "${params.macs_genome_type}" -eq "hs" ]; then
              echo "nopeak	1	15000000	16000000	." >> ${peaks_saf}
          else
              echo "nopeak	chr1	15000000	16000000	." >> ${peaks_saf}
          fi
      fi
    rm "${filt_bam}"
    """

  stub:
    basename = "${experiment_name}_${antibody_name}"
    merged_peaks = "${basename}_peaks_merged.bed"
    peaks_saf = "${basename}_peaks_merged.saf"
    """
    touch "${merged_peaks}" "${peaks_saf}"
    """
}

/*
 * This module defines a process to merge multiple SAF files 
 *
 *
 */

process merge_saf {

  input:
      file(saf_files)
      val(base_name)

  output:
      file(merged_saf_file)

  script:
      merged_saf_file = "${base_name}"+'_merged_peak.saf'
      
      // assume each count file has a header starting with "GeneID" column name
      """
      cat $saf_files | awk 'FNR!=1 && \$1=="GeneID" {next;}{print}' > $merged_saf_file

      """
  stub:
      merged_saf_file = "${base_name}"+'_merged_peak.saf'
      
      """
      touch "${merged_saf_file}"

      """
} 

/*
 * This module defines a process for computing basic ChipSeq QC metrics.
 *
 * Config-defined parameters
 * ----------------------------
 * HOME_REPO - location of the repository
 */
process chip_qc {
  // conda params.HOME_REPO + '/nf/envs/bwa.yaml'

  input:
    tuple val(chipfile_id), file(bam_file), file(saf_file)
    file HOME_REPO

  output:
    tuple val(chipfile_id), file(cell_stats), file(sample_stats)

  script:
    cell_stats = "${chipfile_id}.chipQC_cell.txt"
    sample_stats = "${chipfile_id}.chipQC_sample.txt"

    """
    python "${HOME_REPO}/py/chipQC.py" "${bam_file}" "${saf_file}" "${cell_stats}" --sample_out "${sample_stats}"
    """

  stub:
    cell_stats = "${chipfile_id}.chipQC_cell.txt"
    sample_stats = "${chipfile_id}.chipQC_sample.txt"
    """
    touch "${cell_stats}" "${sample_stats}"
    """
}

/*
 * This module defines a process for merging ChipSeq QC metrics and
 * producing QC plots
 *
 * Config-defined parameters
 * -----------------------------
 * HOME_REPO - location of the repository
 */
process merge_chip_qc {
  // conda params.HOME_REPO + '/nf/envs/skbio.yaml'

  input:
    file chipqc_sample  // expected that .collect() is run. 
    file chipqc_cell // expected that .collect() is run
    val output_base
    file HOME_REPO

  output:
    file chipseq_merged_sample
    file chipseq_merged_cell
    file chipseq_plots

  script:
    chipseq_merged_sample = "${output_base}.sample_chipQC.txt"
    chipseq_merged_cell = "${output_base}.cell_chipQC.txt"
    chipseq_plots = "${output_base}.chipQC.pdf"
    """
    hdr=0
    find . -name '*_cell.txt' > infiles
    while read inf; do
        if [ "\${hdr}" -eq "0" ]; then
            cat \$inf > "${chipseq_merged_cell}"
            hdr=1
        else
            tail -n +2 \$inf >> "${chipseq_merged_cell}"
        fi
    done < infiles
    hdr=0
    find . -name '*_sample.txt' > infiles
    while read inf; do
        if [ "\${hdr}" -eq "0" ]; then
            cat \$inf > "${chipseq_merged_sample}"
            hdr=1
        else
            tail -n +2 \$inf >> "${chipseq_merged_sample}"
        fi
    done < infiles
    python "${HOME_REPO}/py/chipqc_plots.py" "${chipseq_merged_cell}" "${chipseq_merged_sample}" "${chipseq_plots}"
    """

  stub:
    chipseq_merged_sample = "${output_base}.sample_chipQC.txt"
    chipseq_merged_cell = "${output_base}.cell_chipQC.txt"
    chipseq_plots = "${output_base}.chipQC.pdf"
    """
    touch "${chipseq_plots}" "${chipseq_merged_sample}" "${chipseq_merged_cell}"
    """
}
 
/*
 * This module defines a process that, given per-antibody bam files and
 * an output base name, provides a set of peak plots across particular
 * regions of interest.
 * 
 * Config-defined parameters
 * -----------------------------
 * HOME_REPO - location of the repository
 */
process plot_peaks_in_regions {
  conda params.HOME_REPO + "/nf/env/samplot.yaml"

  input:
    file antibody_bams  // .collect() has been run
    file antibody_peaks // .collect() has been run
    val antibody_names // .collect() has been run
    file regions_of_interest
    file genome_gff_gz
    file genome_gff_tbi
    file enhancer_bed_gz
    file enhancer_bed_tbi
    val run_name

  output:
    file './*.png'

  script:
    peak_gz = antibody_peaks.map({ it.toString() + '.gz'})
    """
    # gzip the peak files
    bgzip ${antibody_peaks}
    # tabix them
    for bgp in `ls *.bed.gz`; do
      tabix \$bgp
    done

    # index the bam files
    for bf in `ls *.bam`; do
      samtools index \$bf
    done

    while read line; do
      echo \$line
      chr=\$(echo \$line | awk '{print \$1}')
      str=\$(echo \$line | awk '{print \$2}')
      end=\$(echo \$line | awk '{print \$3}')
      gn=\$(echo \$line | awk '{print \$4}')
    
      samplot plot -n ${antibody_names} \
        -b ${antibody_bams} \
        -c \$chr \
        -s \$str \
        -e \$end \
        -A ${peak_gz} \${enhancer_bed_gz} \
        -T "${genome_gff_gz}" \ 
        -o "${run_name}_\${gn}.png"
    done < regions_of_interest.bed
    """
  
  stub:
    """
    while read line; do
      echo \$line
      gn=\$(echo \$line | awk '{print \$3}')
      touch "${run_name}_\${gn}.png"
    done
    """
}
