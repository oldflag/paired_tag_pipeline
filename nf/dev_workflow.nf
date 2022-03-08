nextflow.enable.dsl=2

/*
 * Development workflow
 */

// CHANGE THIS FILE TO RUN DIFFERENT SAMPLES THROUGH THE PIPELINE
// the digest_file is a csv of the form <sequence_id>,<library_type>,<fastq1>,<fastq2>
// and drives the run of the pipeline

// general parameters
params.RUN_NAME = 'Ishii_ET_RD'
params.HOME_REPO = '/home/app.dev1/repos/pipelines/'
params.py_dir = params.HOME_REPO + 'py/'

// input and output
// LIBRARY_DIGEST = file(params.HOME_REPO + '/ex/example_library_digest.csv')
// SAMPLE_DIGEST = file(params.HOME_REPO + '/ex/zhu2020_sample_digest.csv')
LIBRARY_DIGEST = file('libdigest.csv')
SAMPLE_DIGEST = file('sample_digest.csv')
params.output_dir = 'publisher'


// parameters of R1 trimming
params.trim_ncores = 2
params.adapter_seq = "CTGTCTCTTATA"  // nextera
params.trim_qual = 20

// parameters of R2 parsing
params.linker_file = file(params.HOME_REPO + '/config/linkers.fa')
params.combin_barcodes = file(params.HOME_REPO + '/config/well_barcode_8bp.fa')
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
params.genome_gtf_collapsed_file = file('/home/share/storages/2T/genome/mouse/gencode.vM28.annotation.collapsed.gtf')
params.genome_bed_file = file('/home/share/storages/2T/genome/mouse/gencode.vM28.annotation.bed')
params.genome_element_db = file('/home/chartl/projects/2022-02/annotation_files_for_dna/epimap/mmSDB.saf')
params.count_ncores = 3
params.macs_genome_type = "mm"

// modules
include { trim_fq_single } from params.HOME_REPO + '/nf/modules/trim'
include { parse_pairedtag_r2; split_annot_r1; add_tags as tag1; add_tags as tag2} from params.HOME_REPO + '/nf/modules/pairedtag_reads'
include { star_aligner_single; bwa_aligner_single } from params.HOME_REPO + '/nf/modules/alignment'
include { rnaseqc_call } from params.HOME_REPO + '/nf/modules/rnaseqc'
include { annotate_multiple_features as rna_annot; umitools_count as rna_count; umitools_count as rna_bin_count;
          merge_counts as dna_merge_read; merge_counts as dna_merge_umi; 
          merge_counts as rna_merge_read; merge_counts as rna_merge_umi; 
          merge_counts as rna_merge_bin_read; merge_counts as rna_merge_bin_umi;
          merge_counts as peak_merge_read; merge_counts as peak_merge_umi } from params.HOME_REPO + '/nf/modules/count'
include { annotate_multiple_features as dna_annot; umitools_count as dna_count} from params.HOME_REPO + '/nf/modules/count'
include { annotate_reads_with_features as peak_annot; umitools_count as peak_count } from params.HOME_REPO + '/nf/modules/count'
include { merge_bams as merge_dnabams; merge_bams as merge_annodnabams; merge_bams as merge_annornabams } from params.HOME_REPO + '/nf/modules/alignment'
include { MACS2_peakcall; merge_saf; chip_qc } from params.HOME_REPO + '/nf/modules/peaks'
include { publishData as publishdnabam; publishData as publishrnabam; 
          publishData as publishdnareadcount; publishData as publishdnaumicount; 
          publishData as publishrnareadcount; publishData as publishrnaumicount; 
          publishData as publishrnabinreadcount; publishData as publishrnabinumicount;
          publishData as publishdnapeakreadcount; publishData as publishdnapeakumicount;
          publishData as publishrnaqc } from params.HOME_REPO + '/nf/modules/publish' 

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
  
  //rna_qc
  rnaqc = rnaseqc_call(rna_rawbam.map{it -> tuple(it[0], it[1], it[3], it[4])}, params.genome_gtf_collapsed_file, params.genome_bed_file )

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
  
  /* read and umi count with umi_tools based on a given tag */
  // DNA read and umi count per cell 
  dna_counts = dna_count(dna_withGN[0], 'BN')
  // RNA read and umi count per cell per gene
  rna_counts = rna_count(rna_withGN[0], 'GN')
  // RNA read and umi count per cell  
  rna_bin_counts = rna_bin_count(rna_withGN[0], 'BN')

  
  /* Peak calling with anibodies */
  // grouping by antibody - this is the 5th element of the tuple
  dna_withGN_ab = dna_withGN[0].map{ it -> tuple(it[4], it[1])}.groupTuple()
  // merging bams per antibody
  dna_mg = merge_dnabams(dna_withGN_ab.map{it -> tuple(it[1].collect(), params.RUN_NAME+'_dna_', it[0])})
  // peak calling per antibody and adding antibody name to peak name
  dna_peaks = MACS2_peakcall(dna_mg)

  // bams come out as (bam, expname, antibody)
  // peaks come out as (exp, bed, saf, antibody)
  // desire: (antibody, expname, bam, saf)
  bam_peak_ch = dna_peaks.map{ it -> tuple(it[3], it[0], it[2])}.join(
                   dna_mg.map{it -> tuple(it[2], it[1], it[0])}).map{
                     it -> tuple(it[1] + '_' + it[0], it[4], it[2], 'SAF')}

  chip_qc = chip_qc(bam_peak_ch.map{ it -> tuple(it[0], it[1], it[2])})
  dna_withPeak = peak_annot(bam_peak_ch, 'peaks')


  // combining all peaks for publication
  merged_saf = merge_saf(dna_peaks.map{it -> it[2]}.collect(), "all_antibodies")  
  
  // read and umi count based on peaks
  peak_counts = peak_count(dna_withPeak[0].map{it -> tuple(it[0], it[2], '', '', '')}, 'XT')

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
  read_input = dna_counts[0].map{it -> tuple(1, it[0], it[1], it[2], it[3])}.groupTuple().map{ it -> tuple(it[4].collect(), 'DNA_Q30_aligned_READcount_percell')}
  umi_input = dna_counts[0].map{it -> tuple(1, it[0], it[1], it[2], it[3])}.groupTuple().map{ it -> tuple(it[3].collect(), 'DNA_Q30_aligned_UMIcount_percell')}
  dna_read_merged_h5ad = dna_merge_read(read_input)
  dna_umi_merged_h5ad = dna_merge_umi(umi_input)
  
  
  // merge RNA read and umi counts per cell per gene
  r_read_input = rna_counts[0].map{it -> tuple(1, it[0], it[1], it[2], it[3])}.groupTuple().map{ it -> tuple(it[4].collect(), 'RNA_Q30_aligned_READcount_perGene')}
  r_umi_input = rna_counts[0].map{it -> tuple(1, it[0], it[1], it[2], it[3])}.groupTuple().map{ it -> tuple(it[3].collect(), 'RNA_Q30_aligned_UMIcount_perGene')}
  rna_read_merged_h5ad = rna_merge_read(r_read_input)
  rna_umi_merged_h5ad = rna_merge_umi(r_umi_input)

  // merge RNA read and umi counts per cell
  //TODO: need to group by libary id, assay id and antibody
  r_bin_read_input = rna_bin_counts[0].map{it -> tuple(1, it[0], it[1], it[2], it[3])}.groupTuple().map{ it -> tuple(it[4].collect(), 'RNA_Q30_aligned_READcount_perBIN')}
  r_bin_umi_input = rna_bin_counts[0].map{it -> tuple(1, it[0], it[1], it[2], it[3])}.groupTuple().map{ it -> tuple(it[3].collect(), 'RNA_Q30_aligned_UMIcount_perBIN')}
  rna_bin_read_merged_h5ad = rna_merge_bin_read(r_bin_read_input)
  rna_bin_umi_merged_h5ad = rna_merge_bin_umi(r_bin_umi_input)

  
  // merge Peak read and umi counts
  p_read_input = peak_counts[0].map{ it -> tuple(1, it[0], it[1], it[2], it[3])}.groupTuple().map{ it -> tuple(it[4].collect(), 'DNA_Q30_aligned_READcount_perPeak')}
  p_umi_input = peak_counts[0].map{ it -> tuple(1, it[0], it[1], it[2], it[3])}.groupTuple().map{ it -> tuple(it[3].collect(), 'DNA_Q30_aligned_UMIcount_perPeak')}
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

  //publish QCs 
  publishrnaqc(rnaqc.map{it -> it[3]})
  publishrnabinreadcount(rna_bin_read_merged_h5ad)
  publishrnabinumicount(rna_bin_umi_merged_h5ad)

}
  
