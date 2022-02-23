nextflow.enable.dsl=2

/*
 * Development workflow
 */

// CHANGE THIS FILE TO RUN DIFFERENT SAMPLES THROUGH THE PIPELINE
// the digest_file is a csv of the form <sequence_id>,<library_type>,<fastq1>,<fastq2>
// and drives the run of the pipeline


// general parameters
params.RUN_NAME = 'dev-workflow-run'
params.HOME_REPO = '/home/chartl/repos/pipelines/'
params.py_dir = params.HOME_REPO + 'py/'

// input and output
LIBRARY_DIGEST = file(params.HOME_REPO + '/ex/example_library_digest.csv')
SAMPLE_DIGEST = file(params.HOME_REPO + '/ex/zhu2020_sample_digest.csv')
params.output_dir = '/NAS1/test_runs/2022-02-20/'


// parameters of R1 trimming
params.trim_ncores = 2
params.adapter_seq = "CTGTCTCTTATA"  // nextera
params.trim_qual = 20

// parameters of R2 parsing
params.linker_file = file(params.HOME_REPO + '/config/linkers.fa')
params.combin_barcodes = file(params.HOME_REPO + '/config/well_barcode_7bp.fa')
params.sample_barcodes = SAMPLE_DIGEST  /* now BCs are exported from LIMS file(params.HOME_REPO +  '/config/sample_barcode_4bp.fa') */
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
include { annotate_multiple_features as rna_annot; umitools_count as rna_count; 
          merge_counts as dna_merge_read; merge_counts as dna_merge_umi; 
          merge_counts as rna_merge_read; merge_counts as rna_merge_umi;
          merge_counts as peak_merge_read; merge_counts as peak_merge_umi } from params.HOME_REPO + '/nf/modules/count'
include { annotate_multiple_features as dna_annot; umitools_count as dna_count} from params.HOME_REPO + '/nf/modules/count'
include { annotate_reads_with_features as peak_annot; umitools_count as peak_count } from params.HOME_REPO + '/nf/modules/count'
include { merge_bams as merge_dnabams; merge_bams as merge_annodnabams; merge_bams as merge_annornabams } from params.HOME_REPO + '/nf/modules/alignment'
include { MACS2_peakcall; merge_saf } from params.HOME_REPO + '/nf/modules/peaks'
include { publishData as publishdnabam; publishData as publishrnabam; 
          publishData as publishdnareadcount; publishData as publishdnaumicount; 
          publishData as publishrnareadcount; publishData as publishrnaumicount; 
          publishData as publishdnapeakreadcount; publishData as publishdnapeakumicount } from params.HOME_REPO + '/nf/modules/publish' 

/* channel over rows of the digest */
read1_ch = Channel.fromPath(LIBRARY_DIGEST).splitCsv(header: true, sep: ',')
             .map{ row -> tuple(row.sequence_id, file(row.fastq1)) }
read2_ch = Channel.fromPath(LIBRARY_DIGEST).splitCsv(header: true, sep: ',')
             .map{ row -> tuple(row.sequence_id, file(row.fastq2)) }
type_ch = Channel.fromPath(LIBRARY_DIGEST).splitCsv(header:true, sep: ',')
             .map{ row -> tuple(row.sequence_id, row.library_type) }


workflow {
  
  parsed_barcodes = parse_pairedtag_r2(read2_ch)
  trimmed_reads = trim_fq_single(read1_ch)
  trim_bc_join = trimmed_reads.map{ r -> tuple(r[0], r[1]) }.join(parsed_barcodes)
  run_fastqs = split_annot_r1(trim_bc_join)
   i=0
   j=0
   fastq_ids = run_fastqs[0].map{ it -> [i++, it] }
   fastq_subfiles = run_fastqs[1].map{ it -> [j++, it]}
   // join them
   fqjoin = fastq_ids.join(fastq_subfiles).map{ it -> tuple(it[1], it[2])}
   fqjoin = fqjoin.join(type_ch).map{ it -> tuple(it[0], it[2], it[1])}
   fqjoin = fqjoin.transpose()

   //split into DNA and RNA
   rna_fq = fqjoin.filter{ it[1] =~ /rna/ }.map{it -> tuple(it[0], it[2], it[1])}
   dna_fq = fqjoin.filter{ it[1] =~ /dna/ }.map{it -> tuple(it[0], it[2], it[1])}
  
  //alignment - this will produce a channel of (seq_id, bamfile, seq_type, assay_id, antibody_name)
  dna_rawbam = bwa_aligner_single(dna_fq)
  rna_rawbam = star_aligner_single(rna_fq)
  
  //filtering and tagging
  dna_tagged = tag1(dna_rawbam.filter{ ! it[1].simpleName.contains('_unlinked') })
  rna_tagged = tag2(rna_rawbam.filter{ ! it[1].simpleName.contains('_unlinked') })
  
  // Bin and genome element/gene annotation with featurecounts
  dna_withGN = dna_annot(dna_tagged,
                         tuple(params.genome_bin_file, 'SAF', 'BN'),
                         tuple(params.genome_element_db, 'SAF', 'RE'))
  rna_withGN = rna_annot(rna_tagged,
                         tuple(params.genome_bin_file, 'SAF', 'BN'),
                         tuple(params.genome_gtf_file, 'GTF', 'GN'))
  
  // read and umi count with umi_tools based on a given tag
  dna_counts = dna_count(dna_withGN[0], 'BN')
  rna_counts = rna_count(rna_withGN[0], 'GN')
  
  /* Peak calling with anibodies */
  // grouping by antibody - this is the 5th element of the tuple
  dna_withGN_ab = dna_withGN[0].map{ it -> tuple(it[4], it[1])}.groupTuple()
  // merging bams per antibody
  dna_mg = merge_dnabams(dna_withGN_ab.map{it -> tuple(it[1].collect(), params.RUN_NAME+'_dna_', it[0])})
  // peak calling per antibody and adding antibody name to peak name
  dna_peaks = MACS2_peakcall(dna_mg)

  // combining all peak calling  @hklim  do we want to do this?
  merged_saf = merge_saf(dna_peaks.map{it -> it[2]}.collect(), "all_antibodies")
  // peak annotation
  dna_withPeak = peak_annot(dna_withGN[0].map{it -> tuple(it[0], it[1])},
                            merged_saf,
                            'SAF', 'peaks')
  // read and umi count based on peaks
  peak_counts = peak_count(dna_withPeak[0].map{it -> tuple(it[0], it[2], '')}, 'XT')

  // merge annotated bams
  // note that it[1] here is always 'peaks'
  dna_mg_input = dna_withPeak[0].map{it -> tuple(it[1], it[0], it[1], it[2])}.groupTuple().map{
     it -> tuple(it[3].collect(), params.RUN_NAME + '_anno_', 'dna')
  }
  annodna_mg = merge_annodnabams(dna_mg_input)
  rna_mg_input = rna_withGN[0].map{it -> tuple(it[2], it[0], it[1], it[2])}.groupTuple().map{
     it -> tuple(it[2].collect(), params.RUN_NAME + '_anno_', 'rna')
  }
  annorna_mg = merge_annornabams(rna_mg_input)
  
  // merge DNA read and umi counts
  read_input = dna_counts[0].map{it -> tuple(1, it[0], it[1], it[2], it[3])}.groupTuple().map{ it -> tuple(it[4].collect(), 'DNA_read')}
  umi_input = dna_counts[0].map{it -> tuple(1, it[0], it[1], it[2], it[3])}.groupTuple().map{ it -> tuple(it[3].collect(), 'DNA_umi')}
  dna_read_merged_h5ad = dna_merge_read(read_input)
  dna_umi_merged_h5ad = dna_merge_umi(umi_input)
  
  
  // merge RNA read and umi counts
  r_read_input = rna_counts[0].map{it -> tuple(1, it[0], it[1], it[2], it[3])}.groupTuple().map{ it -> tuple(it[4].collect(), 'RNA_read')}
  r_umi_input = rna_counts[0].map{it -> tuple(1, it[0], it[1], it[2], it[3])}.groupTuple().map{ it -> tuple(it[3].collect(), 'RNA_umi')}
  rna_read_merged_h5ad = rna_merge_read(r_read_input)
  rna_umi_merged_h5ad = rna_merge_umi(r_umi_input)
  
  // merge Peak read and umi counts
  p_read_input = peak_counts[0].map{ it -> tuple(1, it[0], it[1], it[2], it[3])}.groupTuple().map{ it -> tuple(it[4].collect(), 'dna_peak_read')}
  p_umi_input = peak_counts[0].map{ it -> tuple(1, it[0], it[1], it[2], it[3])}.groupTuple().map{ it -> tuple(it[3].collect(), 'dna_peak_umi')}
  peak_read_merged_h5ad = peak_merge_read(p_read_input)
  peak_umi_merged_h5ad = peak_merge_umi(p_umi_input)
  
  // publish results
  publishdnabam(annodna_mg[0].map{ it -> it[0]})
  publishrnabam(annorna_mg[0].map{ it -> it[0]})
  publishdnareadcount(dna_read_merged_h5ad)
  publishdnaumicount(dna_umi_merged_h5ad)
  publishrnareadcount(rna_read_merged_h5ad)
  publishrnaumicount(rna_umi_merged_h5ad)
  publishdnapeakreadcount(peak_read_merged_h5ad)
  publishdnapeakumicount(peak_umi_merged_h5ad)
}
  
