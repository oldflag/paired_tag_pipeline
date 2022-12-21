/*
 * Test workflow for read count merge
 */

nextflow.enable.dsl=2


// parameters
params.HOME_REPO = '/home/app.dev1/repos/pipelines/'
params.datadir = '/home/app.dev1/unittest/mergeCount_h5ad/' 
params.output_dir = '/home/app.dev1/unittest/mergeCount_h5ad/'

umi_count_ch = Channel.fromPath( "${params.datadir}" +'counts/*umiCount.txt.gz')
read_count_ch = Channel.fromPath("${params.datadir}" +'counts/*readCount.txt.gz')

// modules
include { merge_counts } from params.HOME_REPO + '/nf/modules/count'

// Testing module with sample read & umi count files from umi_tools
workflow {

    read_umi_counts = merge_counts(read_count_ch.collect(), umi_count_ch.collect(), "test")
        
}