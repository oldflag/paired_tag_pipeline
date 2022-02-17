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
  conda params.HOME_REPO  + '/nf/envs/macs2.yaml'

  input:
    tuple file(bam_file), val(experiment_name),  val(antibody_name)

  output:
    tuple val(experiment_name), file(merged_peaks), file(peaks_saf)

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
    cat "${merged_peaks}" | awk '{print "${antibody_name}_Peak"NR"\t"\$1"\t"\$2"\t"\$3"\t."}' >> "${peaks_saf}"
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


