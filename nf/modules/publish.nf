/*
 * Modules for publishing data with publishDir
 */
nextflow.enable.dsl=2

/*
 * Publish data to a directory
 * 
 * Config-defined parameters:
 * ----------------------------
 * params.output_dir  - location for publishing data
 *
 */
process publishData {
  
  publishDir "${params.output_dir}", mode: 'copy', overwrite: true

  input:
    file(files)

  output:
    file(files)

  script:
    """ """

  stub:
    """ """
}