/*
 * Modules for parsing / annotating / handling PairedTag reads
 */

/*
 * Downsample fastQ files to the appropriate number of reads
 *
 * Config-defined parameters
 * ---------------------------
 *
 */
process downsample_fastqs {

  input: 
    tuple val(sequence_id), file(r1_fastq), file(r2_fastq), val(library_id)
    val nreads

  output:
    tuple val(sequence_id), file(r1_downsample), file(r2_downsample), val(library_id)

  script:
    r1_downsample = "${r1_fastq}" - ".gz" - ".fq" - ".fastq" + "_ds.fq.gz"
    r2_downsample = "${r2_fastq}" - ".gz" - ".fq" - ".fastq" + "_ds.fq.gz"
    k = nreads * 4

    """
    zcat "${r1_fastq}" | head -n "${k}" | gzip -c > "${r1_downsample}"
    zcat "${r2_fastq}" | head -n "${k}" | gzip -c > "${r2_downsample}"

    """

  stub:
    r1_downsample = "${r1_fastq}" - ".gz" - ".fq" - ".fastq" + "_ds.fq.gz"
    r2_downsample = "${r2_fastq}" - ".gz" - ".fq" - ".fastq" + "_ds.fq.gz"
    """
    touch "${r1_downsample}" "${r2_downsample}"

    """ 
}

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

  // conda params.HOME_REPO + '/nf/envs/skbio.yaml'

  input:
    tuple val(sequence_id), file(r2_fastq), val(library_id), val(library_type)
    


  output:
    tuple val(sequence_id), file(barcode_csv)

  script:
    py_dir = file(params.py_dir)
    combin_barcodes = file(params.combin_barcodes)
    sample_barcodes = file(params.sample_barcodes)
    linker_file = file(params.linker_file)
    barcode_csv = "${sequence_id}.barcodes.csv.gz"

    """
    python "${py_dir}"/parse_R2.py --library_id "${library_id}" --threads "${params.r2_parse_threads}" --umi_size "${params.umi_len}" "${r2_fastq}" "${combin_barcodes}" "${sample_barcodes}" "${linker_file}" "${barcode_csv}"
    """

  stub:
    barcode_csv = "${sequence_id}.barcodes.csv.gz"
    """
    touch "${barcode_csv}"
    """
}


/*
 * Parse R2, annotate R1, and split in a single process.
 *
 * A single sub-library for Paired-Tag (R1.fq R2.fq) will be
 * split into the individual assays (typically 12 tubes, so
 * 12 annotated fq files).
 *
 *
 * * Config-defined parameters
 * ---------------------------
 *  + py_dir: path to the root python directory of the 'pipelines' repo
 *  + combin_barcodes: path to the combinatorial ("well") barcodes sequence file
 *  + sample_barcodes: path to the sample barcodes sequence file - may be a fasta or csv
 *  + linker_file: path to the linker sequence file
 *  + r2_parse_threads: number of threads to use
 *  + umi_len: the UMI length
 */
process process_pairedtag {
  // conda params.HOME_REPO + '/nf/envs/skbio.yaml'

  input:
    tuple val(sequence_id), file(r1_fastq), file(r2_fastq), val(library_id)
    file py_dir
    file combin_barcodes
    file sample_barcodes
    file linker_file

  output:
    val sequence_id
    file 'out/*.fq.gz'
    file sequence_logo

  script:
    sequence_logo="${sequence_id}.logo.pdf"
    
    

    """
    mkdir -p out
    python ${py_dir}/pull_linkers.py "${r2_fastq}" "${sequence_logo}"
    python ${py_dir}/split_pairedtag.py "${r1_fastq}" "${r2_fastq}" "${combin_barcodes}" "${sample_barcodes}" "${linker_file}" ./out/ --library_id "${library_id}" --threads "${params.r2_parse_threads}" --sequence_id "${sequence_id}" --umi_size "${params.umi_len}"
    #for fqf in `ls out`; do
    #    gzip "out/\${fqf}" &
    #done
    #wait 
    """

  stub:
    sequence_logo="${sequence_id}.logo.pdf"
    
    """
    mkdir -p out
    touch "${sequence_logo}"
    touch "out/${sequence_id}__assay01__ab01__sampleX__1.fq.gz"
    touch "out/${sequence_id}__assay02__ab01__sampleY__2.fq.gz"
    touch "out/${sequence_id}__assay01__ab01__sampleX__3.fq.gz"
    touch "out/${sequence_id}__assay02__ab02__sampleY__4.fq.gz"
    """
}


/*
 *
 * Produce a library-level diagnostic plot of barcode parsing
 * 
 * Config-defined parameters
 * ----------------------------
 *  + py_dir: path to the root python directory of the 'pipelines' repo
 *
 */
process barcode_qc {
  // conda params.HOME_REPO + '/nf/envs/skbio.yaml'
  
  input:
    tuple val(sequence_id), file(tagged_fastq)  // using .collect()
    file py_dir
    file plate_layout

  output:
    tuple val(sequence_id), file(barcode_qc_pdf)

  script:
    barcode_qc_pdf = "${sequence_id}.barcode_qc.pdf"

    """
    python ${py_dir}/barcode_qc.py ${tagged_fastq} --output_pdf "${barcode_qc_pdf}" --plate_layout "${plate_layout}"
    """

  stub:
    barcode_qc_pdf = "${sequence_id}.barcode_qc.pdf"
    """
    touch "${barcode_qc_pdf}"
    """
}


