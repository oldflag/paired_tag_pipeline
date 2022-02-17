/*
 * Test workflow for read count merge
 */

nextflow.enable.dsl=2


// parameters
DIGEST = file('run_digest.csv')

// general parameters
params.RUN_NAME = 'dev-workflow-run'
params.HOME_REPO = '/home/app.dev1/repos/pipelines/'
params.py_dir = params.HOME_REPO + 'py/'
params.output_dir = '/home/app.dev1/unittest2/test_parse_r2/'

// parameters of R2 parsing
params.linker_file = file(params.HOME_REPO + '/config/linkers.fa')
params.sample_barcodes = file(params.HOME_REPO +  '/config/sample_barcode_4bp.fa')
params.combin_barcodes = file(params.HOME_REPO + '/config/well_barcode_7bp.fa')
params.umi_len = 10
params.r2_parse_threads = 4


read2_ch = Channel.fromPath(DIGEST).splitCsv(header: true, sep: ',')
             .map{ row -> tuple(row.sequence_id, file(row.read2)) }

// modules
include { parse_pairedtag_r2 } from params.HOME_REPO + '/nf/modules/pairedtag_reads'
include { publishData as publishR2parse} from params.HOME_REPO + '/nf/modules/publish' 

// Testing module with sample read & umi count files from umi_tools
workflow {

    parsed_barcodes = parse_pairedtag_r2(read2_ch)
    publishR2parse(parsed_barcodes.map{ it -> it[1]})
        
}