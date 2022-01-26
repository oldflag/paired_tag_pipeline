/*
 * Modules for RNA and DNA alignment
 */

nextflow.enable.dsl=2

/*
 * This module defines a process for STAR alignment
 *
 * Config-defined parameters:
 * --------------------------
 *   + fastq_trimmed1: read sequence
 *   + fastq_trimmed2: read sequence
 *   + star_index: index
 *   + gtf : annotation file
 *   + alignment_ncore: the number of cores to use for alignment
 */


process star_aligner_single {
  conda params.HOME_REPO + '/envs/star.yaml'
  
    input:
      tuple val(sequence_id), val(fastq_trimmed1)
	  

	output:
	  path ('*.bam'), emit: star_aligned
	  path "*.out", emit: alignment_report
	  file "*SJ.out.tab"
	  file "*Log.out"
	  file "*.bam.bai"
	  file "*.txt"


  script: 
    prefix = "${sequence_id}"
	bamPrefix =  params.outputdir + prefix
	input_fq1 = params.datadir + "${fastq_trimmed1}"

    """
    STAR --readFilesIn $input_fq1 \\
	    --runThreadN $params.alignment_ncore \\
	    --twopassMode Basic \\
	    --outSAMtype BAM SortedByCoordinate \\
	    --readFilesCommand zcat \\
	    --outSAMunmapped Within \\
        --limitBAMsortRAM 16000000000 \\
		--genomeDir $params.star_index \\
	    --outFileNamePrefix $bamPrefix

    """

  stub:
	outbam = bamPrefix + '.Aligned.sortedByCoord.out.bam'
    """
    touch "${outbam}"
    """
}

process star_aligner_pair {
  conda params.HOME_REPO + '/envs/star.yaml'
  
    input:
      tuple val(sequence_id), val(fastq_trimmed1), val(fastq_trimmed2)
	  

	output:
	  path ('*.bam'), emit: star_aligned
	  path "*.out", emit: alignment_report
	  file "*SJ.out.tab"
	  file "*Log.out"
	  file "*.bam.bai"
	  file "*.txt"


  script: 
    prefix = "${sequence_id}"
	bamPrefix =  params.outputdir + prefix
	input_fq1 = params.datadir + "${fastq_trimmed1}"
	input_fq2 = params.datadir + "${fastq_trimmed2}"

    """
    STAR --readFilesIn $input_fq1 $input_fq2 \\
	    --runThreadN "${params.alignment_ncore}" \\
	    --twopassMode Basic \\
	    --outSAMtype BAM SortedByCoordinate \\
	    --readFilesCommand zcat \\
	    --outSAMunmapped Within \\
        --limitBAMsortRAM 16000000000 \\
		--genomeDir $params.star_index \\
	    --outFileNamePrefix $bamPrefix

    """

  stub:
    outbam = bamPrefix + '.Aligned.sortedByCoord.out.bam'
    """
    touch "${outbam}"
    """
}

