nextflow.enable.dsl=2

params.HOME_REPO = '/home/chartl/repos/pipelines/'
params.dataset_name = 'test_sra_file'
params.ncbi_fqdump = '/NAS1/data/paired-tag-2020/sratoolkit.3.0.0-ubuntu64/bin/fasterq-dump'
params.sra_file = file('/home/chartl/projects/2022-05/sra_example.csv')

include { sync_sra } from params.HOME_REPO + '/nf/modules/data'

workflow {
  sync_sra(params.sra_file, params.dataset_name)
}
