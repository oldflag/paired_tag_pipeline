/*
 * Modules for obtaining or otherwise manipulating data
 */
nextflow.enable.dsl=2

/*
 * Given an SRA digest file, synchronize the fastq files
 * to the target directory, update the SRA digest with fastq
 * entries, and return both the finalized digest as well
 * as the individual fastq files
 *
 * Config-defined paramters:
 * ------------------------------
 * ncbi_fqdump - path to the `fasterq-dump` binary
 * output_dir - the publishing directory
 */
process sync_sra {
  input:
    file sra_digest
    val dataset_name

  output:
    file '*.fastq.gz'
    file updated_digest

  script:
    updated_digest = "${dataset_name}_updated.csv"
    """
    IFS="," header=\$(head -n 1 "${sra_digest}")
    i_run=\$(echo \$header | tr ',' '\\n' | nl | grep -w Run | tr -d ' ' | awk -F " " '{print \$1}')
    echo "dataset_name,fastq1,fastq2,fastq3" > extra_cols.txt
    mkdir -p tmp
    while IFS="," read -r run_id; do
      $params.ncbi_fqdump "\${run_id}" -O tmp -e 4 -S -3 2>&1 >> sra.log
      gzip tmp/*
      nfiles=\$(ls -l tmp | grep gz | wc -l)
      files=\$(ls -l tmp | awk '{print \"$params.output_dir/\"\$9}' | grep gz | tr '\\n' ',' | sed 's/,\$//g')
      if [ "\${nfiles}" -eq "1" ]; then
          files="\${files},.,."
      elif [ "\${nfiles}" -eq "2" ]; then
          files="\${files},."
      fi
      echo "${dataset_name},\${files}" >> extra_cols.txt
      mv tmp/*.gz .
    done < <(cut -d "," -f\${i_run} $sra_digest | tail -n +2)
    paste -d , $sra_digest extra_cols.txt > "${dataset_name}_updated.csv"
    """

  stub:
    updated_digest = "${dataset_name}_updated.csv"
    """
    IFS="," header=\$(head -n 1 "${sra_digest}")
    i_run=\$(echo \$header | tr ',' '\\n' | nl | grep -w Run | tr -d ' ' | awk -F " " '{print \$1}')
    echo "dataset_name,fastq1,fastq2,fastq3" > extra_cols.txt
    mkdir -p tmp
    while IFS="," read -r run_id; do
      touch "tmp/\${run_id}_1.fastq.gz" "tmp/\${run_id}_2.fastq.gz" 
      nfiles=\$(ls -l tmp | grep gz | wc -l)
      files=\$(ls -l tmp | awk '{print \"$params.output_dir/\"\$9}' | grep gz | tr '\\n' ',' | sed 's/,\$//g')
      if [ "\${nfiles}" -eq "1" ]; then
          files="\${files},.,."
      elif [ "\${nfiles}" -eq "2" ]; then
          files="\${files},."
      fi
      echo "${dataset_name},\${files}" >> extra_cols.txt
      mv tmp/*.fastq.gz .
    done < <(cut -d "," -f\${i_run} $sra_digest | tail -n +2)
    paste -d , $sra_digest extra_cols.txt > "${dataset_name}_updated.csv"
    """
}
