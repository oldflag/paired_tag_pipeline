# pipelines

Repository for Paired-Tag, RNA-seq, ATAC-seq, ChIP-seq, and other NGS data processing pipelines.

# Pipeline Run with Docker Env
## Local docker container env
### No spikein example
nextflow run ~/repos/pipelines/nf/workflows/prod_pipeline.nf   
-c ~/repos/pipelines/nf/workflows/config/prod_nextflow.config   
--RUN_NAME 'docker_v3_nosplikein_070192023'   
--LIBRARY_DIGEST_FILE '20221108_xxx_MiniSeq_LibraryDigest.csv'   
--SAMPLE_DIGEST_FILE '20221108_xxx_MiniSeq_SampleDigest.csv' 
--SPECIES 'mm' 
--output_dir 'publisher_xxx'  
--NO_SPIKEIN 'yes' 

### Spikein example
nextflow run ~/repos/pipelines/nf/workflows/prod_pipeline.nf   
-c ~/repos/pipelines/nf/workflows/config/prod_nextflow.config   
--RUN_NAME 'dockerV3_spikein_R43_AD_May30' 
--LIBRARY_DIGEST_FILE 'R43_May30.library_digest.csv'   
--SAMPLE_DIGEST_FILE 'R43_May30.sample_digest.csv' 
--NO_SPIKEIN 'no' 
--SPECIES 'hs' 
--SPIKEIN_SPECIES 'mm'    
--output_dir 'publisher_xxx'  
 
