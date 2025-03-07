nextflow.enable.dsl=2

/*
 * PairedTag Pipeline - V2
 */

INSTRUCTIONS = """
//INSTRUCTIONS

//PLEASE COPY AND PASTE THE FOLLOWING INTO THE FILE nextflow.config

executor {
  queueSize=6 
}
conda.enabled = true

conda {
  createTimeout = "120m"
  cacheDir = "/cond"
}

params.RUN_NAME="<your run name>"
params.LIBRARY_DIGEST_FILE="<your library digest>.csv"
params.SAMPLE_DIGEST_FILE="<your sample digest>.csv"
params.SPECIES="mm"  // replace with "hs" for human or "rn" for rat
params.SPIKEIN_SPECIES="hs"  // replase with "mm" for mouse or "rn" for rat
 
params.output_dir="/NAS1/test_runs/" // set to your desired output directory
params.HOME_REPO="/home/chartl/repos/pipelines/" // set to the location of the pipelines repository
"""

def varExists(vne) {
  try {
    vne()
  } catch (exc) {
    false
  }
}

if ( ! varExists({params.SAMPLE_DIGEST_FILE}) ) {
    println INSTRUCTIONS
    System.exit(2)
}

if ( params.RUN_NAME == "<your run name>" ) {
    println "Please remember to fill in the values of run name, library digest, and sample digest"
    System.exit(2)
}


/*** CHANGE NOTHING BELOW HERE ***/

params.py_dir = params.HOME_REPO + "py/"
params.LIBRARY_DIGEST = file(params.LIBRARY_DIGEST_FILE)
params.SAMPLE_DIGEST = file(params.SAMPLE_DIGEST_FILE)

// parameters of R1 trimming
params.trim_ncores = 2
params.adapter_seq = "CTGTCTCTTATA"  // nextera
params.universal_seq = "AGATCGGAAGAG"
params.transposase_seq = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
params.trim_qual = 20

// parameters of R2 parsing
params.linker_file = file(params.HOME_REPO + "/config/linkers.fa")
params.combin_barcodes = file(params.HOME_REPO + "/config/well_barcode_v2_8bp.fa")
params.plate_layout = file(params.HOME_REPO + '/config/plate_layout_v2.csv')
params.sample_barcodes = params.SAMPLE_DIGEST
params.umi_len = 10
params.r2_parse_threads = 2
params.fragment_ncore = 4

params.genome_name = [ 
  "hs": "GRCh38",
  "mm": "GRCm39"
]

params.genome_reference = [
  "hs": file("/home/share/storages/2T/genome/human/GRCh38.primary_assembly.genome.fa"),
  "mm": file("/home/share/storages/2T/genome/mouse/GRCm39.primary_assembly.genome.fa.gz")
]

params.star_index = [
  "hs": file("/home/share/storages/2T/genome/human/star_index/")
  "mm": file("/home/share/storages/2T/genome/mouse/star_index/")
]

params.genome_bin_file = [
  "hs": file("/home/share/storages/2T/genome/human/GRCh38_5kb.saf"),
  "mm": file("/home/share/storages/2T/genome/mouse/GRCm39_5kb.saf")
]

params.genome_saf_file = [
  "hs": file("/home/share/storages/2T/genome/human/gencode.v39.annotation.saf"),
  "mm": file("/home/share/storages/2T/genome/mouse/gencode.vM28.annotation.saf")
]

params.genome_gtf_collapsed_file = [
  "hs": file("/home/share/storages/2T/genome/human/gencode.v39.annotation.collapsed.gtf"),
  "mm": file("/home/share/storages/3T/genome/mouse/gencode.vM28.annotation.collapsed.gtf")
]

params.heterochromatin_saf_file = [
  "hs": file("/home/chartl/projects/2022-08/GBM/GSE127123_GBM_heterochromatin.hg39.saf"),
  "mm": file("/home/share/storages/2T/genome/mouse/20220603-mm39-brain-heterochromatin.saf")
]

params.promoter_saf_file = [
  "hs": file("/home/share/storages/2T/genome/human/RegEmtDB_promoter_hg38.saf"),
  "mm": file("/home/share/storages/2T/genome/mouse/old_annots/GRCm39_Encode_Promoters.saf")
]

params.enhancer_saf_file = [
  "hs": file("/home/share/storages/2T/genome/human/RegEmtDB_enhancer_hg38.saf"),
  "mm": file("/home/share/storages/2T/genome/mouse/old_annots/GRCm39_Encode_Enhancers.saf")
]

params.alignment_ncore = 4
params.fragment_ncore = 6
params.ramsize = 5000000000
params.count_ncores = 3
params.fragment_ncores = 4


// modules
include { trim_fq_single as trim_rna; trim_fq_single as trim_dna} from params.HOME_REPO + "/nf/modules/trim"
include { process_pairedtag; 
          barcode_qc as barcode_qc;
          add_tags as tag_rna; 
          add_tags as tag_dna;
          add_tags as tag_rna_spike;
          add_tags as tag_dna_spike } from params.HOME_REPO + "/nf/modules/pairedtag_reads"
include { star_aligner_single; bwa_aligner_single } from params.HOME_REPO + "/nf/modules/alignment_spikein"
include { alignment_qc as alignment_qc_spike; 
          alignment_qc; 
          merge_alignment_qc as merge_aligment_qc_spike;
          merge_alignment_qc
          bam_to_frag } from params.HOME_REPO + "/nf/modules/alignment"
include { rnaseqc_call; merge_rnaseqc } from params.HOME_REPO + "/nf/modules/rnaseqc_spikein"
include { umitools_count as rna_count; umitools_count as rna_bin_count;
          umitools_count as dna_count;
          merge_counts as dna_merge_read; merge_counts as dna_merge_umi; 
          merge_counts as rna_merge_read; merge_counts as rna_merge_umi; 
          merge_counts as rna_merge_bin_read; merge_counts as rna_merge_bin_umi;
          merge_counts as peak_merge_read; merge_counts as peak_merge_umi} from params.HOME_REPO + "/nf/modules/count"
include { annotate_reads_with_features as peak_annot; umitools_count as peak_count } from params.HOME_REPO + "/nf/modules/count"
include { h5ad_qc; cluster_qc } from params.HOME_REPO + "/nf/modules/count"
include { merge_bams as merge_dna_bams; merge_bams as merge_rna_bams } from params.HOME_REPO + "/nf/modules/alignment"
include { MACS2_multi; MACS2_multi as MACS2_multi_spike; 
          merge_saf; merge_saf as merge_saf_spike; 
          chip_qc; chip_qc as chip_qc_spike; 
          merge_chip_qc; merge_chip_qc as merge_chip_qc_spike } from params.HOME_REPO + "/nf/modules/peaks"
include { publishData as publishdnabam; publishData as publishrnabam; 
          publishData as publishdnareadcount; publishData as publishdnaumicount; 
          publishData as publishrnareadcount; publishData as publishrnaumicount; 
          publishData as publishrnabinreadcount; publishData as publishrnabinumicount;
          publishData as publishdnapeakreadcount; publishData as publishdnapeakumicount;
          publishData as publishfragments; catAndPublish as publishlogs;
          publishData as publishlogo;
          publishData as publishrnaqc;
          publishData as publishbarcodeqc;
          publishData as publishalignmentqc;
          publishData as publishchipqc;
          publishData as publishh5adqc;
          publishData as publishclusterqc;
          publishData as publish10xrna;
          publishData as publish10xdna;
          publishData as publishtracks } from params.HOME_REPO + "/nf/modules/publish" 

/* channel over rows of the digest */
pair_ch = Channel.fromPath(params.LIBRARY_DIGEST).splitCsv(header: true, sep: ",").map{ row -> tuple(row.sequence_id, file(row.fastq1), file(row.fastq2), row.lysis_id, row.library_type)}
type_ch = Channel.fromPath(params.LIBRARY_DIGEST).splitCsv(header:true, sep: ",").map{ row -> tuple(row.sequence_id, row.library_type) }

def inner_join(ch_a, ch_b) {
    return ch_b.cross(ch_a).map { [it[0][0], *it[1][1..-1], *it[0][1..-1]] }
}

workflow {
  // # note: params.SAMPLE_DIGEST is used implicitly here
  //pair_ch.subscribe{ println(it) }
  split_fqs = process_pairedtag(pair_ch.map{ it -> tuple(it[0], it[1], it[2], it[3])})
  //split_fqs[0].subscribe{ println(it) }
  publishlogo(split_fqs[2])
  i=0
  j=0
  fastq_ids = split_fqs[0].map{ it -> [i++, it] }
  fastq_subfiles = split_fqs[1].map{ it -> [j++, it]}
   // join them
  fqjoin = fastq_ids.join(fastq_subfiles).map{ it -> tuple(it[1], it[2])}
  fqjoin = fqjoin.join(type_ch).map{ it -> tuple(it[0], it[2], it[1])}
  fqjoin = fqjoin.transpose()
  // fqjoin.subscribe{println(it)}
   // this is now a tuple of (seq_id, seq_type, fastq)
  //fqjoin.map{it -> tuple(it[0], it[2])}.groupTuple().subscribe{println it}
  barcode_pdfs = barcode_qc(fqjoin.map{ it -> tuple(it[0], it[2]) }.groupTuple())

  publishbarcodeqc(barcode_pdfs)
   
  rna_fq = fqjoin.filter{ it[1] =~ /rna/ }.map{it -> tuple(it[0], it[2], it[1])}
  dna_fq = fqjoin.filter{ it[1] =~ /dna/ }.map{it -> tuple(it[0], it[2], it[1])}

  // in v2 pipeline, a single .bam is produced
  // in v3 pipeline, a (filtered) primary bam is produced in slot 0
  //                 a (filtered) spike-in bam i sproduced in slot 1
  //                 a spike-in info file is produced in slot 2

  dna_raw_bams = bwa_aligner_single(dna_fq)
  rna_raw_bams = star_aligner_single(rna_fq)

  primary_bams_all = dna_raw_bams[0].mix(rna_raw_bams[0])
  spikein_bams_all = dna_raw_bams[1].mix(rna_raw_bams[1])

  bamqc = alignment_qc(all_bams)
  bamqc_spike = alignment_qc_spike(spikein_bams_all)
  alignment_qcfile = merge_alignment_qc(bamqc.map{ it -> it[1]}.collect(), params.RUN_NAME)
  alignment_qcfile_spike = merge_alignment_qc(bamqc_spike.map{ it -> it[1]}.collect(), params.RUN_NAME + '_spikein')
  publishalignmentqc(alignment_qcfile.mix(alignment_qcfile_spike))

  // fragments
  fragfiles = bam_to_frag(dna_rawbam[0].mix(dna_rawbam[1]).map{ it -> it[1] })  // [0]: fragments [1]: index [2]: logs
  publishfragments(fragfiles[0])
  publishlogs(fragfiles[2].collect(), "gather_convert_fragments.log")


  //rna_qc
  rnaqc = rnaseqc_call(rna_raw_bams[0].mix(rna_raw_bams[1]).map{it -> tuple(it[0], it[1], it[3], it[4])}, params.genome_gtf_collapsed_file)
  rnaqc_mg = merge_rnaseqc(rnaqc.map{it -> it[4]}.collect(), params.RUN_NAME)

  
  // Peak calling with antibodies //
  // grouping by antibody - this is the 5th element of the tuple
  dna_byAB = dna_raw_bam[0].map{ it -> tuple(it[4], it[1])}.groupTuple().map{ it -> tuple(it[0], it[1], params.RUN_NAME)}
  dna_byAB_spike = dna_raw_bam[1].map{ it -> tuple(it[4], it[1])}.groupTuple().map{ it -> tuple(it[0], it[1], params.RUN_NAME)}
  
  // peak calling per antibody and adding antibody name to peak name
  
  // combining all peaks for publication
  dna_peaks = MACS2_multi(dna_byAB)  // input: (antibody_name, bam_file_list, experiment_name)
  dna_peaks_spike = MACS2_multi(dna_byAB_spike)  // input: (antibody_name, bam_file_list, experiment_name)

  // output: (experiment, bed, saf, antibody)
  // output[1]: wig file
  merged_saf = merge_saf(dna_peaks[0].map{it -> it[2]}.collect(), "all_antibodies")
  merged_saf_spike = merge_saf_spike(dna_peaks_spike[0].map{ it -> it[2] }.collect(), "_spike_all_antibodies")
  publishtracks(dna_peaks[1].mix(dna_peaks_spike[1]).collect())


  //filtering and tagging
  dna_tagged = tag_dna(dna_raw_bam[0].filter{ ! it[1].simpleName.contains("_unlinked") },
                       tuple(params.genome_bin_file[params.SPECIES], "BN", "SAF"),
                       tuple(params.enhancer_saf_file[params.SPECIES], "EE", "SAF"),
                       tuple(params.promoter_saf_file[params.SPECIES], "EP", "SAF"),
                       tuple(params.genome_saf_file[params.SPECIES], "GN", "SAF"),
                       merged_saf.map{it -> tuple(it, "PK", "SAF")},
                       tuple(file("input.6"), "NA", "NAN"))

  dna_tagged_spike = tag_dna_spike(dna_raw_bam[1].filter{ ! it[1].simpleName.contains("_unlinked")},
                       tuple(params.genome_bin_file[params.SPIKEIN_SPECIES], "BN", "SAF"),
                       tuple(params.enhancer_saf_file[params.SPIKEIN_SPECIES], "EE", "SAF"),
                       tuple(params.promoter_saf_file[params.SPIKEIN_SPECIES], "EP", "SAF"),
                       tuple(params.genome_saf_file[params.SPIKEIN_SPECIES], "GN", "SAF"),
                       merged_saf_spike.map{it -> tuple(it, "PK", "SAF")},
                       tuple(file("input.6"), "NA", "NAN"))


  // val(alignment_id), file(annot_bam_file), val(seqtype), val(assay_id), val(antibody)
  // (chipfile_id, bam_file, saf_file)
  // create (antibody, assay, sample, bam, antibody_peak)
  // strategy:  (antibody, assay, alignment, bam).join(antibody, saf).map(assay"_"antibody, alignment, bam, saf)
  chipqc_input = inner_join(dna_tagged.map{ it -> tuple(it[4], it[3], it[0], it[1])},
    dna_peaks[0].map{ it -> tuple(it[3], it[2]) }
  ).map{it -> tuple(it[2] + "_" + it[0] + "_" + it[1], it[3], it[4])} // use the alignment id
  //chipqc_input.subscribe{ println(it) }
  chip_qc = chip_qc(chipqc_input)
  chip_qc_mg = merge_chip_qc(chip_qc.map{it -> it[2]}.collect(),
                             chip_qc.map{it -> it[1]}.collect(), params.RUN_NAME)

  chipqc_input_spike = inner_join(dna_tagged_spike.map{ it -> tuple(it[4], it[3], it[0], it[1])},
    dna_peaks_spike[0].map{ it -> tuple(it[3], it[2]) }
  ).map{it -> tuple(it[2] + "_" + it[0] + "_" + it[1], it[3], it[4])} // use the alignment id
  //chipqc_input.subscribe{ println(it) }
  chip_qc_spike = chip_qc_spike(chipqc_input_spike)
  chip_qc_mg_spike = merge_chip_qc_spike(chip_qc.map{it -> it[2]}.collect(),
                             chip_qc.map{it -> it[1]}.collect(), params.RUN_NAME)

  rna_tagged = tag_rna(rna_rawbam.filter{ ! it[1].simpleName.contains("_unlinked") },
                       tuple(params.genome_bin_file[params.SPECIES], "BN", "SAF"),
                       tuple(params.genome_saf_file[params.SPECIES], "GN", "SAF"),
                       tuple(file("input.3"), "NA", "NAN"),
                       tuple(file("input.4"), "NA", "NAN"),
                       tuple(file("input.5"), "NA", "NAN"),
                       tuple(file("input.6"), "NA", "NAN"))
 
  
  // read and umi count with umi_tools based on a given tag //
  // DNA read and umi count per cell 
  dna_counts = dna_count(dna_tagged, "BN")
  // RNA read and umi count per cell per gene
  rna_counts = rna_count(rna_tagged, "GN")
  // RNA read and umi count per cell  
  rna_bin_counts = rna_bin_count(rna_tagged, "BN")
  
  // read and umi count based on peaks
  peak_counts = peak_count(dna_tagged, "PK")


  
  // merge DNA read and umi counts
  read_input = dna_counts[0].map{it -> tuple(1, it[0], it[1], it[2], it[3])}.groupTuple().map{ it -> tuple(it[4].collect(), "DNA_Q30_aligned_READcount_percell")}
  umi_input = dna_counts[0].map{it -> tuple(1, it[0], it[1], it[2], it[3])}.groupTuple().map{ it -> tuple(it[3].collect(), "DNA_Q30_aligned_UMIcount_percell")}
  dna_read_merged_h5ad = dna_merge_read(read_input, params.SAMPLE_DIGEST)
  dna_umi_merged_h5ad = dna_merge_umi(umi_input, params.SAMPLE_DIGEST)
  
  
  // merge RNA read and umi counts per cell per gene
  r_read_input = rna_counts[0].map{it -> tuple(1, it[0], it[1], it[2], it[3])}.groupTuple().map{ it -> tuple(it[4].collect(), "RNA_Q30_aligned_READcount_perGene")}
  r_umi_input = rna_counts[0].map{it -> tuple(1, it[0], it[1], it[2], it[3])}.groupTuple().map{ it -> tuple(it[3].collect(), "RNA_Q30_aligned_UMIcount_perGene")}
  rna_read_merged_h5ad = rna_merge_read(r_read_input, params.SAMPLE_DIGEST)
  rna_umi_merged_h5ad = rna_merge_umi(r_umi_input, params.SAMPLE_DIGEST)

  // merge RNA read and umi counts per cell
  //TODO: need to group by libary id, assay id and antibody
  r_bin_read_input = rna_bin_counts[0].map{it -> tuple(1, it[0], it[1], it[2], it[3])}.groupTuple().map{ it -> tuple(it[4].collect(), "RNA_Q30_aligned_READcount_perBIN")}
  r_bin_umi_input = rna_bin_counts[0].map{it -> tuple(1, it[0], it[1], it[2], it[3])}.groupTuple().map{ it -> tuple(it[3].collect(), "RNA_Q30_aligned_UMIcount_perBIN")}
  rna_bin_read_merged_h5ad = rna_merge_bin_read(r_bin_read_input, params.SAMPLE_DIGEST)
  rna_bin_umi_merged_h5ad = rna_merge_bin_umi(r_bin_umi_input, params.SAMPLE_DIGEST)

  
  // merge Peak read and umi counts
  p_read_input = peak_counts[0].map{ it -> tuple(1, it[0], it[1], it[2], it[3])}.groupTuple().map{ it -> tuple(it[4].collect(), "DNA_Q30_aligned_READcount_perPeak")}
  p_umi_input = peak_counts[0].map{ it -> tuple(1, it[0], it[1], it[2], it[3])}.groupTuple().map{ it -> tuple(it[3].collect(), "DNA_Q30_aligned_UMIcount_perPeak")}
  peak_read_merged_h5ad = peak_merge_read(p_read_input, params.SAMPLE_DIGEST)
  peak_umi_merged_h5ad = peak_merge_umi(p_umi_input, params.SAMPLE_DIGEST)

  // qc the final h5ad
  h5qc_out = h5ad_qc(rna_read_merged_h5ad[0], rna_umi_merged_h5ad[0], dna_read_merged_h5ad[0], dna_umi_merged_h5ad[0])
  h5_pdf = h5qc_out[0]

  clustqc_out = cluster_qc(h5qc_out[1], h5qc_out[2])

  // convert scanpy-conforming h5ad files into 10X/Seurat-conforming h5 files
  //rna_seur = convert_rna_umi(rna_umi_merged_h5ad[0], params.genome_name)
  //dna_seur = convert_dna_umi(dna_umi_merged_h5ad[0], params.genome_name)
  
  // publish results
  publishdnabam(dna_tagged.map{ it -> it[1]})
  publishrnabam(rna_tagged.map{ it -> it[1]})
  publishdnareadcount(dna_read_merged_h5ad[0])
  publishdnaumicount(dna_umi_merged_h5ad[0])
  publishrnareadcount(rna_read_merged_h5ad[0])
  publishrnaumicount(rna_umi_merged_h5ad[0])
  publishdnapeakreadcount(peak_read_merged_h5ad[0])
  publishdnapeakumicount(peak_umi_merged_h5ad[0])

  //publish QCs 
  publishrnaqc(rnaqc_mg[1])
  publishrnabinreadcount(rna_bin_read_merged_h5ad[0])
  publishrnabinumicount(rna_bin_umi_merged_h5ad[0])
  publishchipqc(chip_qc_mg[2])
  publishh5adqc(h5_pdf)
  publishclusterqc(clustqc_out)
  //publish10xrna(rna_seur)
  //publish10xdna(dna_seur)
}
  
