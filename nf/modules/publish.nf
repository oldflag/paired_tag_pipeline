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


/*
 * Publish files to a directory, but first concatenate them
 *
 * Config-defined parameters:
 * -------------------------------
 * params.output_dir - location for publishing data
 */
process catAndPublish {
   publishDir "${params.output_dir}", mode: 'copy', overwrite: true

   input:
     file(files)
     val(out_fname)

   output:
     file outfile

   script:
     outfile=out_fname
     """
     touch catAndPublishRun
     cat $files > $out_fname
     """

   stub:
     outfile=out_fname
     """
     touch $out_fname
     """
}    

/*
 * Concatenate multiple dataframes and publish
 */
process catAndPublishDF {
  publishDir "${params.output_dir}", mode: 'copy', overwrite: true

  input:
    file(files)
    val(out_fname)

  output:
    file outfile

  script:
    outfile=out_fname
    """
    touch catAndPublishRun
    python "${params.HOME_REPO}/py/cat_df.py" $files --out $outfile
    """

  stub:
    outfile=out_fname
    """
    touch $out_fname
    """
}
