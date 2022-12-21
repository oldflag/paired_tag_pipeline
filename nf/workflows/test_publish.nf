/*
 * Test workflow for read count merge
 */

nextflow.enable.dsl=2


// parameters
params.HOME_REPO = '/home/app.dev1/repos/pipelines/'
params.datadir = '/home/app.dev1/tmp_h5ad_merge/' 
params.output_dir = '/home/app.dev1/unittest2/publisher/'

umi_count_ch = Channel.fromPath( "${params.datadir}" +'counts/*umiCount.txt.gz')
read_count_ch = Channel.fromPath("${params.datadir}" +'counts/*readCount.txt.gz')

// modules
include { publishData } from params.HOME_REPO + '/nf/modules/publish'

// Testing module with sample bam and count files
workflow {

    publishData(read_count_ch.collect())
   // publishData(umi_count_ch.collect())
        
}