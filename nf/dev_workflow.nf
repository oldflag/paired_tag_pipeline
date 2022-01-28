nextflow.enable.dsl=2

/*
 * Development workflow
 */

// CHANGE THIS FILE TO RUN DIFFERENT SAMPLES THROUGH THE PIPELINE
// the digest_file is a csv of the form <sequence_id>,<read1_fq>,<read2_fq>,key1=value1;key2=value2;<etc>
// and drives the run of the pipeline
DIGEST = file('run_digest.csv')

// general parameters
params.HOME_REPO = '/home/chartl/repos/pipelines/'
params.py_dir = params.HOME_REPO + 'py/'
params.output_dir = '/home/chartl/projects/2022-01/nextflow-dev/'

// parameters of R1 trimming
params.trim_ncores = 2
params.adapter_seq = "CTGTCTCTTATA"  // nextera
params.trim_qual = 20

// parameters of R2 parsing
params.linker_file = file(params.HOME_REPO + '/config/linkers.fa')
params.sample_barcodes = file(params.HOME_REPO +  '/config/sample_barcode_4bp.fa')
params.combin_barcodes = file(params.HOME_REPO + '/config/well_barcode_7bp.fa')
params.umi_len = 10
params.r2_parse_threads = 4

// modules
include { trim_fq_single } from params.HOME_REPO + '/nf/modules/trim'
include { parse_pairedtag_r2; split_annot_r1 } from params.HOME_REPO + '/nf/modules/pairedtag_reads'


/* channel over rows of the digest */
read1_ch = Channel.fromPath(DIGEST).splitCsv(header: true, sep: ',')
             .map{ row -> tuple(row.sequence_id, row.read1) }
read2_ch = Channel.fromPath(DIGEST).splitCsv(header: true, sep: ',')
             .map{ row -> tuple(row.sequence_id, row.read2) }

workflow {
  parsed_barcodes = parse_pairedtag_r2(read2_ch)
  trimmed_reads = trim_fq_single(read1_ch)

  parsed_barcodes.subscribe{ println it}
  trimmed_reads.subscribe{ println it }

  trim_bc_join = trimmed_reads.map{ r -> tuple(r[0], r[1]) }.join(parsed_barcodes)

  trim_bc_join.subscribe{ println it.flatten() }

  run_fastqs = split_annot_r1(trim_bc_join)
}
  
