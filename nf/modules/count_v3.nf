/*
 * Modules for counting reads in genomic locations
 */

nextflow.enable.dsl=2

/*
 * This module defines a process to count UMI within features
 * in a feature-annotated bam file, where the 'XT' feature is
 * used to determine read --> feature assignments. These tags
 * are added by featureCounts.
 *
 * Config-defined parameters:
 * -------------------------
 *   + HOME_REPO : the path to the home repository
 *
 */
process umitools_count {
  // conda params.HOME_REPO + '/nf/envs/umi_tools.yaml'
  input:
      tuple val(sequence_id), file(annot_bam), val(seqtype), val(assay_id), val(antibody)
      val count_tag

  output:
      tuple val(sequence_id), val(count_tag), file(umi_counts_txt), file(read_counts_txt)
      file umi_log

  script:
      bamfn = annot_bam.simpleName
      umi_counts_txt = bamfn - '.bam' + '_' + count_tag + '_umiCount.txt.gz'
      read_counts_txt = bamfn - '.bam' + '_' + count_tag + '_readCount.txt.gz'
      umi_log = bamfn - '.bam' + '_' + count_tag + '_umitools.log'
      """
      samtools index "${annot_bam}"

      umi_tools count \
        --extract-umi-method=tag \
        --umi-tag=MI \
        --cell-tag=CB \
        --method unique \
        --per-gene \
        --gene-tag="${count_tag}" \
        --per-cell \
        --mapping-quality 0 \
        --cell-tag-split=@ \
        -I "${annot_bam}" \
        -S "${umi_counts_txt}" 2>&1 > "${umi_log}"

      umi_tools count \
        --extract-umi-method=tag \
        --umi-tag=XX \
        --cell-tag=CB \
        --method unique \
        --per-gene \
        --gene-tag="${count_tag}" \
        --per-cell \
        --cell-tag-split=@ \
        --mapping-quality 0 \
        -I "${annot_bam}" \
        -S "${read_counts_txt}" 2>&1 >> "${umi_log}"
        
      """

    stub:
      bamfn = annot_bam.simpleName
      umi_counts_txt = bamfn - '.bam' + '_' + count_tag + '_umiCount.txt.gz'
      read_counts_txt = bamfn - '.bam' + '_' + count_tag + '_readCount.txt.gz'
      umi_log = bamfn - '.bam' + '_' + count_tag + '_umitools.log'
      """
      touch "${umi_counts_txt}"
      touch "${read_counts_txt}"
      touch "${umi_log}"
      """
}


/*
 * This module defines a process to merge multiple UMI and Read count files into one respectively 
 * and convert the merged files into H5AD format 
 *
 * Config-defined parameters:
 * -------------------------
 *   + HOME_REPO : the path to the home repository
 *
 */

process merge_counts {
  
  // conda params.HOME_REPO + '/nf/envs/scanalysis.yaml'

  input:
      tuple file(count_files), val(file_header)
      file sample_digest_csv

  output:
      file merged_count_h5ad
      file merged_count

  script:
      merged_count_h5ad = "${file_header}"+'_merged.h5ad'
      merged_count = "${file_header}"+'_merged.txt.gz'
      
      // assume each count file has a header starting with "gene" column name
      """
      zcat $count_files | awk 'FNR!=1 && \$1=="gene" {next;}{print}' | gzip -c > "${merged_count}"

      python "${params.HOME_REPO}/py/count2h5ad.py" "${merged_count}" "${sample_digest_csv}" "${merged_count_h5ad}"


      """
  stub:
      merged_count_h5ad = "${file_header}"+'_merged.h5ad'
      merged_count = "${file_header}"+'_merged.txt.gz'
      """
      touch "${merged_count_h5ad}"
      touch "${merged_count}"
      """
}

/*
 * This module defines a process to produce basic QC outputs from final read / umi counts
 * 
 * Config-defined parameters:
 * -------------------------
 *   + HOME_REPO : the path to the home repository
 *   + RUN_NAME : the name of the run
 */
process h5ad_qc {
  // conda params.HOME_REPO + '/nf/envs/scanalysis.yaml'

  input:
    file rna_reads
    file rna_umi
    file dna_reads
    file dna_umi

  output:
    file out_pdf
    file out_rna
    file out_dna

  script:
    out_pdf = "${params.RUN_NAME}_h5adqc.pdf"
    out_rna = (rna_umi.toString() - '.h5ad' + '.analysis.h5ad')
    out_dna = (dna_umi.toString() - '.h5ad' + '.analysis.h5ad')
    """
    python "${params.HOME_REPO}/py/h5ad_qc.py" "${rna_reads}" "${rna_umi}" "${dna_reads}" "${dna_umi}" "${out_pdf}"
    """

  stub:
    out_pdf = "${params.RUN_NAME}_h5adqc.pdf"
    out_rna = (rna_umi.toString() - '.h5ad' + '.analysis.h5ad')
    out_dna = (dna_umi.toString() - '.h5ad' + '.analysis.h5ad')
    """
    touch $out_pdf $out_rna $out_dna
    """
}


/*
 * This process calls a python script which performs cell selection
 * and clustering on RNA and DNA components, and produces a
 * QC pdf file
 */
process cluster_qc {
  // conda params.HOME_REPO + '/nf/envs/scanalysis.yaml'

  input:
    file rna_umi
    file dna_umi

  output:
    file out_pdf

  script:
    out_pdf = "${params.RUN_NAME}_cluster_qc.pdf"
    """
    python "${params.HOME_REPO}/py/cluster_qc.py" "${dna_umi}" "${rna_umi}" "${out_pdf}"
    """

  stub:
    out_pdf = "${params.RUN_NAME}_cluster_qc.pdf"
    """
    touch $out_pdf
    """
}
