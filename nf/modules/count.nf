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
  conda params.HOME_REPO + '/nf/envs/featurecounts.yaml'
  input:
    tuple val(sequence_id), file(bam_file)
    file annotation_file
    val annotation_type  // 'GTF' or 'SAF'
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

    samtools sort -n "fc_out/${bamname}.featureCounts.bam" > tmp.bam
    samtools sort tmp.bam > "${annot_bam}"
    rm "fc_out/${bamname}.featureCounts.bam" tmp.bam
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
  conda params.HOME_REPO + '/nf/envs/umitools.yaml'
  input:
      tuple val(sequence_id), file(annot_bam)
      val count_tag

  output:
      tuple val(sequence_id), val(count_tag), file(umi_counts), file(read_counts)
      file umi_log

  script:
      bamfn = annot_bam.simpleName
      umi_counts_txt = bamfn - '.bam' + '_' + count_tag + '_umiCount.txt.gz'
      read_counts_txt = bamfn - '.bam' + '_' + count_tag + '_readCount.txt.gz'
      umi_counts = bamfn - '.bam' + '_' + count_tag + '_umiCount.h5ad'
      read_counts = bamfn - '.bam' + '_' + count_tag + '_readCount.h5ad'
      umi_log = bamfn - '.bam' + '_' + count_tag + '_umitools.log'
      """
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

      python "${params.HOME_REPO}/py/count2h5ad.py" "${umi_counts_txt}" "${umi_counts}"

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

      python "${params.HOME_REPO}/py/count2h5ad.py" "${read_counts_txt}" "${read_counts}"
      """

    stub:
      bamfn = annot_bam.simpleName
      umi_counts_txt = bamfn - '.bam' + '_' + count_tag + '_umiCount.txt.gz'
      read_counts_txt = bamfn - '.bam' + '_' + count_tag + '_readCount.txt.gz'
      umi_counts = bamfn - '.bam' + '_' + count_tag + '_umiCount.h5ad'
      read_counts = bamfn - '.bam' + '_' + count_tag + '_readCount.h5ad'
      umi_log = bamfn - '.bam' + '_' + count_tag + '_umitools.log'
      """
      touch "${umi_counts}"
      touch "${read_counts}"
      touch "${umi_log}"
      """
}


process annotate_multiple_features {
  conda params.HOME_REPO + '/nf/envs/featurecounts.yaml'
  input:
    tuple val(sequence_id), file(bam_file)
    tuple file(annotation_file1), val(annotation_type1), val(destination_tag1)
    tuple file(annotation_file2), val(annotation_type2), val(destination_tag2)

  output:
    tuple val(sequence_id), file(merged_bam)
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

    samtools sort -n "fc_out/${bamname}.featureCounts.bam" > tmp.bam
    samtools sort tmp.bam > "${annot_bam1}"
    rm "fc_out/${bamname}.featureCounts.bam" tmp.bam

    featureCounts -F "${annotation_type2}" \
      -O -Q 30 \
      -T "${params.count_ncores}" \
      --verbose -R BAM \
      --Rpath fc_out \
      -a "${annotation_file2}" \
      -o "${count_file2}" \
      "${bam_file}" 2>&1 >> "${fc_log}"

    samtools sort -n "fc_out/${bamname}.featureCounts.bam" > tmp.bam 
    samtools sort tmp.bam > "${annot_bam2}"
    rm "fc_out/${bamname}.featureCounts.bam" tmp.bam

    python "${params.HOME_REPO}/py/combine_tags.py" \ 
      --drop XN,XS \
      "${annot_bam1}:XT:${destination_tag1}"  \
      "${annot_bam2}:XT:${destination_tag2}" \
      "${merged_bam}" 2>&1 > "${merge_log}"
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

