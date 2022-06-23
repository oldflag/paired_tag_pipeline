/*
 * Modules for counting reads in genomic locations
 */

nextflow.enable.dsl=2

/*
 * This module defines a process for annotating a bam file with
 * genomic features using featureCounts from `subread`
 *
 * Config-defined parameters:
 * --------------------------
 *   + HOME_REPO : the path to the home repository
 *   + count_ncores : the number of cores to use for counting
 */
process annotate_reads_with_features {
  // conda params.HOME_REPO + '/nf/envs/featurecounts.yaml'
  input:
    tuple val(sequence_id), file(bam_file), file(annotation_file), val(annotation_type)  // annotation_type : 'SAF' or 'GTF'
    val annotation_name  // e.g., 'genes', 'enhancers', 'bins'

  output:
    tuple val(sequence_id), val(annotation_name), file(annot_bam)
    file fc_log

  script:
    bamname = bam_file.simpleName
    count_file = bamname - '.bam' + '.fc_counts.txt'
    annot_bam = bamname - '.bam' + '_' + annotation_name + '_annot.bam'
    fc_log = bamname - '.bam' + '_' + annotation_name + '_fc.log'

    """
    mkdir -p fc_out
    featureCounts -F "${annotation_type}" \
      -O -Q 30 \
      -T "${params.count_ncores}" \
      --verbose -R BAM \
      --Rpath fc_out \
      -a "${annotation_file}" \
      -o "${count_file}" \
      "${bam_file}" 2>&1 > "${fc_log}"

    samtools sort -n "fc_out/${bamname}.bam.featureCounts.bam" > tmp.bam
    samtools sort tmp.bam > "${annot_bam}"
    rm "fc_out/${bamname}.bam.featureCounts.bam" tmp.bam
    """

  stub:
    bamname = bam_file.simpleName
    count_file = bamname - '.bam' + '.fc_counts.txt'
    annot_bam = bamname - '.bam' + '_annot.bam'
    fc_log = bamname - '.bam' + '_fc.log'

    """
    touch "${annot_bam}"
    touch "${fc_log}"
    """
}


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
        --mapping-quality 30 \
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
        --mapping-quality 30 \
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


process annotate_multiple_features {
  // conda params.HOME_REPO + '/nf/envs/featurecounts.yaml'
  input:
    tuple val(sequence_id), file(bam_file), val(seqtype), val(assay_id), val(antibody)
    tuple file(annotation_file1), val(annotation_type1), val(destination_tag1)
    tuple file(annotation_file2), val(annotation_type2), val(destination_tag2)
    file(HOME_REPO)

  output:
    tuple val(sequence_id), file(merged_bam), val(seqtype), val(assay_id), val(antibody)
    tuple file(fc_log), file(merge_log)

  script:
    bamname = bam_file.simpleName
    count_file1 = bamname - '.bam' + '.fc_count1.txt'
    annot_bam1 = bamname - '.bam' + '_annot1.bam'
    count_file2 = bamname - '.bam' + '.fc_count2.txt'
    annot_bam2 = bamname - '.bam' + '_annot2.bam'
    fc_log = bamname - '.bam' + '_multi_annot.log'
    merge_log = bamname - '.bam' + '_tag_merge.log'
    merged_bam = bamname - '.bam' + '_multiAnnot.bam'
    """
    mkdir -p fc_out
    featureCounts -F "${annotation_type1}" \
      -O -Q 30 \
      -T "${params.count_ncores}" \
      --verbose -R BAM \
      --Rpath fc_out \
      -a "${annotation_file1}" \
      -o "${count_file1}" \
      "${bam_file}" 2>&1 > "${fc_log}"

    samtools sort -n "fc_out/${bamname}.bam.featureCounts.bam" > tmp.bam
    samtools sort tmp.bam > "${annot_bam1}"
    rm "fc_out/${bamname}.bam.featureCounts.bam" tmp.bam

    featureCounts -F "${annotation_type2}" \
      -O -Q 30 \
      -T "${params.count_ncores}" \
      --verbose -R BAM \
      --Rpath fc_out \
      -a "${annotation_file2}" \
      -o "${count_file2}" \
      "${bam_file}" 2>&1 >> "${fc_log}"

    samtools sort -n "fc_out/${bamname}.bam.featureCounts.bam" > tmp.bam 
    samtools sort tmp.bam > "${annot_bam2}"
    rm "fc_out/${bamname}.bam.featureCounts.bam" tmp.bam

    pyf="${HOME_REPO}/py/combine_tags.py"

    samtools index "${annot_bam1}"
    samtools index "${annot_bam2}"

    python "\${pyf}" "${annot_bam1}:XT:${destination_tag1}"  "${annot_bam2}:XT:${destination_tag2}" "${merged_bam}" 2>&1 > "${merge_log}"


    """

  stub:
    bamname = bam_file.simpleName
    count_file1 = bamname - '.bam' + '.fc_count1.txt'
    annot_bam1 = bamname - '.bam' + '_annot1.bam'
    count_file2 = bamname - '.bam' + '.fc_count2.txt'
    annot_bam2 = bamname - '.bam' + '_annot2.bam'
    fc_log = bamname - '.bam' + '_multi_annot.log'
    merge_log = bamname - '.bam' + '_tag_merge.log'
    merged_bam = bamname - '.bam' + '_multiAnnot.bam'

    """
    touch "${merged_bam}"
    touch "${fc_log}"
    touch "${merge_log}"
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
      file HOME_REPO

  output:
      file(merged_count_h5ad)
      file(merged_count)

  script:
      merged_count_h5ad = "${file_header}"+'_merged.h5ad'
      merged_count = "${file_header}"+'_merged.txt.gz'
      
      // assume each count file has a header starting with "gene" column name
      """
      zcat $count_files | awk 'FNR!=1 && \$1=="gene" {next;}{print}' | gzip -c > "${merged_count}"

      python "${HOME_REPO}/py/count2h5ad.py" "${merged_count}" "${sample_digest_csv}" "${merged_count_h5ad}"


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
 *
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
    file HOME_REPO

  output:
    file out_pdf

  script:
    out_pdf = "${params.RUN_NAME}_h5adqc.pdf"
    """
    python "${HOME_REPO}/py/h5ad_qc.py" "${rna_reads}" "${rna_umi}" "${dna_reads}" "${dna_umi}" "${out_pdf}"
    """

  stub:
    out_pdf = "${params.RUN_NAME}_h5adqc.pdf"
    """
    touch $out_pdf
    """
}
