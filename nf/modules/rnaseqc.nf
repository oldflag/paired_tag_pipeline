/*
 * Modules for RNA QC
 */
nextflow.enable.dsl=2


/*
 * This defines a process for merging RNA-SeQC csv outputs
 * into a single file, and generating a pdf
 *
 * Config-defined parameters:
 * ----------------------------
 * HOME_REPO - location of the repository
 */
process merge_rnaseqc {
  conda params.HOME_REPO + '/nf/envs/skbio.yaml'  // piggyback on skbio for plotting packages

  input:
    file metrics_csv  // called with .collect()
    val base_name

  output:
    file merged_metrics
    file metrics_plots

   script:
     merged_metrics = base_name + '_RNASeQC_merged.csv'
     metrics_plots = base_name + '_RNASeQC_merged.pdf'
     """
     find . -name '*.csv' > qclist.txt
     hdr=0
     while read qcf; do
         if [ \$hdr -eq "0" ]; then
             cat "\${qcf}" > "${merged_metrics}"
             hdr=1
         else
             tail -n 1 "\${qcf}" >> "${merged_metrics}"
         fi
     done < qclist.txt

     python "${params.HOME_REPO}/py/plot_rnaseqc.py" "${merged_metrics}" "${metrics_plots}"
     """

   stub:
     merged_metrics = base_name + '_RNASeQC_merged.csv'
     metrics_plots = base_name + '_RNASeQC_merged.pdf'
     """
     touch "${merged_metrics}" "${metrics_plots}"
     """
}

/*
 * This defines a process for RNA QC using RNA-SeQC package
 *
 * Config-defined parameters:
 * ----------------------------
 * HOME_REPO - location of the repository
 * 
 *
 */
process rnaseqc_call {
  conda params.HOME_REPO  + '/nf/envs/rseqc.yaml'

  input:
    tuple val(sequence_id), file(bam_file), val(assay), val(antibody_name) 
    file(gtf_file)

  output:
    tuple val(sequence_id), val(assay), val(antibody_name), file(metrics_tsv), file(metrics_csv), file(exon_reads_gct), file(gene_reads_gct), file(gene_tpm_gct), file(gene_fragments_gct), file(coverage_tsv)

  script:
    
    basename = bam_file.simpleName
    metrics_tsv = "${basename}.metrics.tsv"
    metrics_csv = "${basename}.metrics.csv"
    exon_reads_gct = "${basename}.exon_reads.gct" 
    gene_reads_gct = "${basename}.gene_reads.gct"
    gene_tpm_gct = "${basename}.gene_tpm.gct" 
    gene_fragments_gct = "${basename}.gene_fragments.gct"
    coverage_tsv = "${basename}.coverage.tsv"


    """
    rnaseqc $gtf_file $bam_file --sample $basename --coverage . --unpaired
    cut -f1 "${metrics_tsv}" | sed 's/, /_/g' | tr '\n' ',' | sed 's/,\$/\\n/g' > "${metrics_csv}"
    cut -f2 "${metrics_tsv}" | sed 's/, /_/g' | tr '\n' ',' | sed 's/,\$/\\n/g' >> "${metrics_csv}"
    """

  stub:

    basename = bam_file.simpleName
    metrics_tsv = "${basename}.metrics.tsv"
    metrics_csv = "${basename}.metrics.csv"
    exon_reads_gct = "${basename}.exon_reads.gct" 
    gene_reads_gct = "${basename}.gene_reads.gct"
    gene_tpm_gct = "${basename}.gene_tpm.gct" 
    gene_fragments_gct = "${basename}.gene_fragments.gct"
    coverage_tsv = "${basename}.coverage.tsv"

    """
    touch ${basename}.metrics.tsv
    touch ${basename}.metrics.csv
    touch ${basename}.exon_reads.gct 
    touch ${basename}.gene_reads.gct
    touch ${basename}.gene_tpm.gct 
    touch ${basename}.gene_fragments.gct
    touch ${basename}.coverage.tsv
    """
}


