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
    tuple val(sequence_id), file(r2_fastq), val(library_id)
    file py_dir
    file combin_barcodes
    file sample_barcodes
    file linker_file

  output:
    tuple val(sequence_id), file(barcode_csv)

  script:
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
    file 'out/*1.fq.gz'
    file 'out/*2.fq.gz'
    file sequence_logo

  script:
    sequence_logo="${sequence_id}.logo.pdf"
    """
    mkdir -p out
    python "${py_dir}/pull_linkers.py" "${r2_fastq}" "${sequence_logo}"
    python "${py_dir}/split_pairedtag_v3.py" "${r1_fastq}" "${r2_fastq}" "${combin_barcodes}" "${sample_barcodes}" "${linker_file}" ./out/ --library_id "${library_id}" --threads "${params.r2_parse_threads}" --sequence_id "${sequence_id}" --umi_size "${params.umi_len}"
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
    touch "out/${sequence_id}__assay01__ab01__sampleX1__1.fq.gz"
    touch "out/${sequence_id}__assay01__ab01__sampleX1__2.fq.gz"
    touch "out/${sequence_id}__assay02__ab01__sampleY1__1.fq.gz"
    touch "out/${sequence_id}__assay02__ab01__sampleY1__2.fq.gz"
    touch "out/${sequence_id}__assay01__ab02__sampleX2__1.fq.gz"
    touch "out/${sequence_id}__assay01__ab02__sampleX2__2.fq.gz"
    touch "out/${sequence_id}__assay02__ab02__sampleY2__1.fq.gz"
    touch "out/${sequence_id}__assay02__ab02__sampleY2__2.fq.gz"
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
    python "${py_dir}"/barcode_qc.py $tagged_fastq --output_pdf "${barcode_qc_pdf}" --plate_layout "${plate_layout}"
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
    file py_dir
    file combin_barcodes
    file sample_barcodes

  output:
  
    val sequence_id
    file 'out/*.fq.gz'

  script:
    """
    mkdir -p out
    python "${py_dir}/annotate_split_R1.py" "${r1_trim_fq}" "${barcode_csv}" "${combin_barcodes}" "${sample_barcodes}" --outdir out --sequence_id "${sequence_id}" --noestimate

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
    annot_bam_file = bam_file.simpleName - '.bam' + '_tag.bam'
    """
    touch $annot_bam_file
    """
}
