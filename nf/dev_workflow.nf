nextflow.enable.dsl=2

/*
 * Development workflow
 */

// CHANGE THIS FILE TO RUN DIFFERENT SAMPLES THROUGH THE PIPELINE
// the digest_file is a csv of the form <sequence_id>,<read1_fq>,<read2_fq>,key1=value1;key2=value2;<etc>
// and drives the run of the pipeline
DIGEST = file('run_digest.csv')

// general parameters
params.RUN_NAME = 'dev-workflow-run'
params.HOME_REPO = '/home/app.dev1/repos/pipelines/'
params.py_dir = params.HOME_REPO + 'py/'
params.output_dir = '/home/app.dev1/unittest2/result/'

// parameters of R1 trimming
params.trim_ncores = 2
params.adapter_seq = "CTGTCTCTTATA"  // nextera
params.trim_qual = 20

// parameters of R2 parsing
params.linker_file = file(params.HOME_REPO + '/config/linkers.fa')
params.sample_barcodes = file(params.HOME_REPO +  '/config/sample_barcode_4bp.fa')
params.combin_barcodes = file(params.HOME_REPO + '/config/well_barcode_7bp.fa')
params.umi_len = 10
params.r2_parse_threads = 4

// parameters of genome alignment
params.genome_reference = file('/home/share/storages/2T/genome/mouse/GRCm39.primary_assembly.genome.fa.gz')
params.alignment_ncore = 4
params.ramsize = 2000000000
params.star_index = file('/home/share/storages/2T/genome/mouse/star_index/')

// parameters for read counting & bam annotation
params.genome_bin_file = file('/home/share/storages/2T/genome/mouse/GRCm39_1kb.saf')
params.genome_gtf_file = file('/home/share/storages/2T/genome/mouse/gencode.vM28.annotation.gtf')
params.genome_element_db = file('/home/chartl/projects/2022-02/annotation_files_for_dna/epimap/mmSDB.saf')
params.count_ncores = 3
params.macs_genome_type = "mm"

// modules
include { trim_fq_single } from params.HOME_REPO + '/nf/modules/trim'
include { parse_pairedtag_r2; split_annot_r1; add_tags as tag1; add_tags as tag2} from params.HOME_REPO + '/nf/modules/pairedtag_reads'
include { star_aligner_single; bwa_aligner_single } from params.HOME_REPO + '/nf/modules/alignment'
include { annotate_multiple_features as rna_annot; umitools_count as rna_count; merge_counts as dna_merge; merge_counts as rna_merge; merge_counts as peak_merge } from params.HOME_REPO + '/nf/modules/count'
include { annotate_multiple_features as dna_annot; umitools_count as dna_count} from params.HOME_REPO + '/nf/modules/count'
include { annotate_reads_with_features as peak_annot; umitools_count as peak_count } from params.HOME_REPO + '/nf/modules/count'
include { merge_bams as merge_dnabams; merge_bams as merge_annodnabams; merge_bams as merge_annornabams } from params.HOME_REPO + '/nf/modules/alignment'
include { MACS2_peakcall } from params.HOME_REPO + '/nf/modules/peaks'
include { publishData as publishdnabam; publishData as publishrnabam; 
          publishData as publishdnareadcount; publishData as publishdnaumicount; 
          publishData as publishrnareadcount; publishData as publishrnaumicount; 
          publishData as publishdnapeakreadcount; publishData as publishdnapeakumicount } from params.HOME_REPO + '/nf/modules/publish' 

/* channel over rows of the digest */
read1_ch = Channel.fromPath(DIGEST).splitCsv(header: true, sep: ',')
             .map{ row -> tuple(row.sequence_id, file(row.read1)) }
read2_ch = Channel.fromPath(DIGEST).splitCsv(header: true, sep: ',')
             .map{ row -> tuple(row.sequence_id, file(row.read2)) }
kv_ch = Channel.fromPath(DIGEST).splitCsv(header:true, sep: ',')
             .map{ row -> tuple(row.sequence_id, row.keyvalues) }

workflow {

  parsed_barcodes = parse_pairedtag_r2(read2_ch)
  trimmed_reads = trim_fq_single(read1_ch)

  trim_bc_join = trimmed_reads.map{ r -> tuple(r[0], r[1]) }.join(parsed_barcodes)

  run_fastqs = split_annot_r1(trim_bc_join)
  // we have to join the two outputs into a single one
  // add indexes
  i=0
  j=0
  fastq_ids = run_fastqs[0].map{ it -> [i++, it] }
  fastq_subfiles = run_fastqs[1].map{ it -> [j++, it]}
  // join them
  fqjoin = fastq_ids.join(fastq_subfiles).map{ it -> tuple(it[1], it[2])}
  fqjoin = fqjoin.join(kv_ch).map{ it -> tuple(it[0], it[2], it[1])}
  fqjoin = fqjoin.transpose()
  rna_fq = fqjoin.filter{ it[1] == 'type=rna' }.map{it -> tuple(it[0], it[2])}
  dna_fq = fqjoin.filter{ it[1] == 'type=dna'}.map{it -> tuple(it[0], it[2])}
  dna_rawbam = bwa_aligner_single(dna_fq)
  rna_rawbam = star_aligner_single(rna_fq)
  dna_tagged = tag1(dna_rawbam.filter{ ! it[1].simpleName.contains('trimmed_unlinked') })
  rna_tagged = tag2(rna_rawbam.filter{ ! it[1].simpleName.contains('trimmed_unlinked') })
  dna_withGN = dna_annot(dna_tagged,
                         tuple(params.genome_bin_file, 'SAF', 'BN'),
                         tuple(params.genome_element_db, 'SAF', 'RE'))
  rna_withGN = rna_annot(rna_tagged,
                         tuple(params.genome_bin_file, 'SAF', 'BN'),
                         tuple(params.genome_gtf_file, 'GTF', 'GN'))

  dna_counts = dna_count(dna_withGN[0], 'BN')
  rna_counts = rna_count(rna_withGN[0], 'GN')
  dna_mg = merge_dnabams(dna_withGN[0].map{ it -> it[1] }.collect(),
                      params.RUN_NAME + '_dna')
  dna_peaks = MACS2_peakcall(dna_mg)

  dna_withPeak = peak_annot(dna_withGN[0],
                            dna_peaks.map{it -> it[2]}.collect(),
                            'SAF', 'peaks')
  peak_counts = peak_count(dna_withPeak[0].map{it -> tuple(it[0], it[2])}, 'XT')
  
  // merge annotated bams

  annodna_mg = merge_annodnabams(dna_withPeak[0].map{ it -> it[2] }.collect(),
                      params.RUN_NAME + '_anno_dna') 
  
  annorna_mg = merge_annornabams(rna_withGN[0].map{ it -> it[1] }.collect(),
                      params.RUN_NAME + '_anno_rna') 
  // merge count results
  dna_merged_h5ad = dna_merge(dna_counts[0].map{ it -> it[3] }.collect(),
                              dna_counts[0].map{ it -> it[2] }.collect(),
                              "DNA")
 
  rna_merged_h5ad = rna_merge(rna_counts[0].map{ it -> it[3] }.collect(),
                              rna_counts[0].map{ it -> it[2] }.collect(),
                              "RNA")


  peak_merged_h5ad = peak_merge(peak_counts[0].map{ it -> it[3] }.collect(),
                                peak_counts[0].map{ it -> it[2] }.collect(),
                                "Peak")


  // publishing

  publishdnabam(annodna_mg[0].map{ it -> it[1]})
  publishrnabam(annorna_mg[0].map{ it -> it[1]})

  publishdnareadcount(dna_merged_h5ad[0])
  publishdnaumicount(dna_merged_h5ad[1])
  publishrnareadcount(rna_merged_h5ad[0])
  publishrnaumicount(rna_merged_h5ad[1])
  publishdnapeakreadcount(peak_merged_h5ad[0])
  publishdnapeakumicount(peak_merged_h5ad[1])

}
  
