/*
 * Modules for parsing / annotating / handling PairedTag reads
 */


/*
 * Apply the paired-tag R2 parser to an R2 of a paired-tag
 * sequence file, producing a CSV
 *
 * Config-defined parameters
 * ---------------------------
 *  + py_dir: path to the root python directory of the 'pipelines' repo
 *  + combin_barcodes: path to the combinatorial ("well") barcodes sequence file
 *  + sample_barcodes: path to the sample barcodes sequence file - may be a fasta or csv
 *  + linker_file: path to the linker sequence file
 *  + r2_parse_threads: number of threads to use
 *  + umi_len: the UMI length
 */
process parse_pairedtag_r2 {

  conda params.HOME_REPO + '/nf/envs/skbio.yaml'

  input:
    tuple val(sequence_id), file(r2_fastq)

  output:
    tuple val(sequence_id), file(barcode_csv)

  script:
    barcode_csv = "${sequence_id}.barcodes.csv.gz"

    """
    python "${params.py_dir}"/parse_R2.py --threads "${params.r2_parse_threads}" --umi_size "${params.umi_len}" "${r2_fastq}" "${params.combin_barcodes}" "${params.sample_barcodes}" "${params.linker_file}" "${barcode_csv}"
    """

  stub:
    barcode_csv = "${sequence_id}.barcodes.csv.gz"
    """
    touch "${barcode_csv}"
    """
}


/*
 *
 * Produce a library-level diagnostic plot of barcode parsing
 * 
 * Config-defined parameters
 * ----------------------------
 *  + py_dir: path to the root python directory of the 'pipeliens' repo
 *
 */
process barcode_qc {
  conda params.HOME_REPO + '/nf/envs/skbio.yaml'
  
  input:
    tuple val(sequence_id), file(barcode_csv)

  output:
    tuple val(sequence_id), file(barcode_qc_pdf)

  script:
    barcode_qc_pdf = "${sequence_id}.barcode_qc.pdf"
    """
    python "${params.py_dir}"/barcode_qc.py "${barcode_csv}" "${barcode_qc_pdf}"
    """

  stub:
    barcode_qc_pdf = "${sequence_id}.barcode_qc.pdf"
    """
    touch "${barcode_qc_pdf}"
    """
}

/*
 * Combine parsed sequence barcodes with trimmed fastq reads, and split these
 * reads into fastQs with ~10K cells in them
 *
 * TODO: Optional arguments for [expected # cells] and [target # cells output]
 *
 * Config-defined parameters
 * ---------------------------
 *  + py_dir: path to the root python directory of the 'pipelines' repo
 *  + combin_barcodes: path to the combinatorial ("well") barcodes sequence file
 *  + sample_barcodes: path to the sample barcodes sequence file
 *  + linker_file: path to the linker sequence file
 */ 
process split_annot_r1 {
  input:
    tuple val(sequence_id), file(r1_trim_fq), file(barcode_csv)

  output:
    val sequence_id
    file 'out/*.fq.gz'

  script:
    """
    mkdir -p out
    python "${params.py_dir}/annotate_split_R1.py" "${r1_trim_fq}" "${barcode_csv}" "${params.combin_barcodes}" "${params.sample_barcodes}" --outdir out --sequence_id "${sequence_id}" --noestimate

    nul=\$(zcat out/*__unknown__unlinked__1.fq.gz | wc -l)
    if [ "\${nul}" -lt "12" ]; then
        rm -f out/*__unknown__unlinked__1.fq.gz
    fi

    nul=\$(zcat out/*__sample__unlinked__1.fq.gz | wc -l)
    if [ "\${nul}" -lt "12" ]; then
       rm -f out/*__sample__unlinked__1.fq.gz
    fi
    """

  stub:
    """
    mkdir -p out
    touch "out/${sequence_id}__assay01__ab01__1.fq.gz"
    touch "out/${sequence_id}__assay02__ab01__2.fq.gz"
    touch "out/${sequence_id}__assay01__ab01__3.fq.gz"
    touch "out/${sequence_id}__assay02__ab02__4.fq.gz"
    """
}

process add_tags {
  conda params.HOME_REPO + '/nf/envs/bwa.yaml'
  
  input:
    tuple val(alignment_id), file(bam_file), val(seqtype), val(assay_id), val(antibody)

  output:
    tuple val(alignment_id), file(annot_bam_file), val(seqtype), val(assay_id), val(antibody)

  script:
    unsorted_bam = bam_file.simpleName - '.bam' + '_tag_uns.bam'
    annot_bam_file = bam_file.simpleName - '.bam' + '_tag.bam'
    """
    python "${params.py_dir}"/add_tags.py "${bam_file}" "${unsorted_bam}" --library "${alignment_id}" --antibody "${antibody}"
    samtools sort "${unsorted_bam}" -o "${annot_bam_file}"
    """

  stub:
    annot_bam_file = bam_file.simpleName.replace('.bam', '_tag.bam')
    """
    touch "${annot_bam_file}"
    """
}
