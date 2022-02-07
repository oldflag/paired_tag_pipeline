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
    tuple val(experiment_name), file(bam_file)

  output:
    tuple val(experiment_name), file(merged_peaks), file(peaks_saf)

  script:
    narrow_peaks = "${experiment_name}_peaks.narrowPeak"
    broad_peaks = "${experiment_name}_peaks.broadPeak"
    merged_peaks = "${experiment_name}_peaks_merged.bed"
    peaks_saf = "${experiment_name}_peaks_merged.saf"
    filt_bam = "${experiment_name}_filt.bam"
    """
    samtools view -hb -q 30 "${bam_file}" > "${filt_bam}" 
    macs2 callpeak -t "${filt_bam}" -n "${experiment_name}" --outdir . \
       -q 0.1 -g "${params.macs_genome_type}" --nomodel
    macs2 callpeak -t "${filt_bam}" --broad -n "${experiment_name}" --outdir . \
       -q 0.1 -g "${params.macs_genome_type}" --nomodel

    cat "${narrow_peaks}" "${broad_peaks}" | cut -f1-5 | bedtools sort -i /dev/stdin | bedtools merge -i /dev/stdin > "${merged_peaks}"
    echo "GeneID	Chr	Start	End	Strand" > "${peaks_saf}"
    cat "${merged_peaks}" | awk '{print "Peak"NR"\t"\$1"\t"\$2"\t"\$3"\t."}' >> "${peaks_saf}"
    rm "${filt_bam}"
    """

  stub:
    merged_peaks = "${experiment_name}_peaks_merged.bed"
    peaks_saf = "${experiment_name}_peaks_merged.saf"
    """
    touch "${merged_peaks}" "${peaks_saf}"
    """
}
    


