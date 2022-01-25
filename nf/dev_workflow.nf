nextflow.enable.dsl=2

/*
 * Development workflow
 */

// CHANGE THIS FILE TO RUN DIFFERENT SAMPLES THROUGH THE PIPELINE
// the digest_file is a csv of the form <sample_id>,<read1_fq>,<read2_fq>,key1=value1;key2=value2;<etc>
// and drives the run of the pipeline
DIGEST = file('run_digest.csv')

// general parameters
HOME_REPO = '/home/chartl/repos/pipelines/'
params.py_dir = HOME_REPO + 'py/'

// parameters of R1 trimming
params.trim_ncores = 2
params.adapter_seq = "CTGTCTCTTATA"  // nextera
params.trim_qual = 20

// parameters of R2 parsing
params.linker_file = file('./config/linkers.fa')
params.sample_barcodes = file('./config/sample_barcode_4bp.fa')
params.combin_barcodes = file('./config/well_barcode_7bp.fa')
params.umi_len = 10

// modules
include { trim_fq_single } from './nf/modules/trim'


/* channel over rows of the digest */
read1_ch = Channel.fromPath(DIGEST).splitCsv(header: true, sep: ',')
             .map{ row -> (row.sample_id, row.read1) }
read2_ch = Channel.fromPath(DIGEST).splitCsv(header: true, sep: ',')
             .map{ row -> (row.sample_id, row.read2) }
