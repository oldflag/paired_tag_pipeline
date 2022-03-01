/*
 * Modules for RNA QC
 */
nextflow.enable.dsl=2

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
    file(bed_file)

  output:
    tuple val(sequence_id), val(assay), val(antibody_name), file(metrics_tsv), file(exon_reads_gct), file(gene_reads_gct), file(gene_tpm_gct), file(gene_fragments_gct), file(coverage_tsv)

  script:
    
    basename = bam_file.simpleName
    metrics_tsv = "${basename}.metrics.tsv"
    exon_reads_gct = "${basename}.exon_reads.gct" 
    gene_reads_gct = "${basename}.gene_reads.gct"
    gene_tpm_gct = "${basename}.gene_tpm.gct" 
    gene_fragments_gct = "${basename}.gene_fragments.gct"
    coverage_tsv = "${basename}.coverage.tsv"


    """
    rnaseqc $gtf_file $bam_file --bed $bed_file --sample $basename --coverage .
    """

  stub:

    basename = bam_file.simpleName
    metrics_tsv = "${basename}.metrics.tsv"
    exon_reads_gct = "${basename}.exon_reads.gct" 
    gene_reads_gct = "${basename}.gene_reads.gct"
    gene_tpm_gct = "${basename}.gene_tpm.gct" 
    gene_fragments_gct = "${basename}.gene_fragments.gct"
    coverage_tsv = "${basename}.coverage.tsv"

    """
    touch ${basename}.metrics.tsv
    touch ${basename}.exon_reads.gct 
    touch ${basename}.gene_reads.gct
    touch ${basename}.gene_tpm.gct 
    touch ${basename}.gene_fragments.gct
    touch ${basename}.coverage.tsv
    """
}


