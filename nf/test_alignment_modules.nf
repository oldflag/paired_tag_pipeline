/*
 * Test workflow for RNA and DNA alignment
 */

nextflow.enable.dsl=2


// parameters for alignment
params.HOME_REPO = '/home/app.dev1/test'
params.datadir = '/home/hklim/projects/tmp/' 
params.star_index = '/home/share/storages/2T/index/human/STAR' 
params.alignment_ncore = 8
params.outputdir = './result/'

test_ch_single = Channel
    .from( 'sequence_id,fastq_trimmed1\nIPmyc_205_control_1,IPmyc_205_control_1.1.fastq.gz.1' )
    .splitCsv(header: true, sep: ',')
	.map{ row -> tuple(row.sequence_id, row.fastq_trimmed1)}

test_ch_pair = Channel
    .from( 'sequence_id,fastq_trimmed1,fastq_trimmed2\nIPmyc_205_control_1,IPmyc_205_control_1.1.fastq.gz.1,IPmyc_205_control_1.2.fastq.gz.1' )
    .splitCsv(header: true, sep: ',')
	.map{ row -> tuple(row.sequence_id, row.fastq_trimmed1, row.fastq_trimmed2)}

// modules
include { star_aligner_single; star_aligner_pair } from params.HOME_REPO + '/modules/alignment'

// Testing modules with single or paired reads
workflow {
  
// star_result = star_aligner_pair(test_ch_pair)  // 
  star_result = star_aligner_single(test_ch_single)
  
}
    
