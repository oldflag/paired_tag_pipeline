# pipelines

Repository for Paired-Tag, RNA-seq, ATAC-seq, ChIP-seq, and other NGS data processing pipelines.

Pipeline Run with Docker Env
example:  
nextflow run /where/is/pipelines/nf/workflows/dev_docker_v2.nf \
-c /where/is/pipelines/nf/workflows/config/dev_v2_nextflow.config \
--RUN_NAME 'test_altos' \
--LIBRARY_DIGEST_FILE '20221108_Altos_MiniSeq_LibraryDigest2.csv' \
--SAMPLE_DIGEST_FILE '20221108_Altos_MiniSeq_SampleDigest2.csv' \
--SPECIES 'mm'
