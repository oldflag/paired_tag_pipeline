nextflow.enable.dsl=2

/*
 * PairedTag Pipeline - V3
 */

INSTRUCTIONS = """
//INSTRUCTIONS
Please see the following configuration file and update as needed.
./config/dev_docker_v3_nextflow.config for DEV
./config/prod_nextflow.config for PROD

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

params.py_dir = file(params.HOME_REPO + 'py')
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
params.sample_barcodes = params.SAMPLE_DIGEST
params.umi_len = 10
params.r2_parse_threads = 2
params.fragment_ncore = 4

params.genome_name = [ 
  "hs": "GRCh38",
  "mm": "GRCm39"
]

params.genome_reference = [
  "hs": file(params.GENOME_DIR + "/human/GRCh38.primary_assembly.genome.fa"),
  "mm": file(params.GENOME_DIR + "/mouse/GRCm39.primary_assembly.genome.fa.gz")
]

params.star_index = [
  "hs": file(params.GENOME_DIR + "/human/star_index/"),
  "mm": file(params.GENOME_DIR + "/mouse/star_index/")
]

params.genome_bin_file = [
  "hs": file(params.GENOME_DIR + "/human/GRCh38_5kb.saf"),
  "mm": file(params.GENOME_DIR + "/mouse/GRCm39_5kb.saf")
]

params.genome_saf_file = [
  "hs": file(params.GENOME_DIR + "/human/gencode.v39.annotation.saf"),
  "mm": file(params.GENOME_DIR + "/mouse/gencode.vM28.annotation.saf")
]

params.genome_gtf_collapsed_file = [
  "hs": file(params.GENOME_DIR + "/human/gencode.v39.annotation.collapsed.gtf"),
  "mm": file(params.GENOME_DIR + "/mouse/gencode.vM28.annotation.collapsed.gtf")
]

params.heterochromatin_saf_file = [
  "hs": file(params.GENOME_DIR + "/human/GSE127123_GBM_heterochromatin.hg39.saf"),
  "mm": file(params.GENOME_DIR + "/mouse/20220603-mm39-brain-heterochromatin.saf")
]

params.promoter_saf_file = [
  "hs": file(params.GENOME_DIR + "/human/RegEmtDB_promoter_hg38.saf"),
  "mm": file(params.GENOME_DIR + "/mouse/old_annots/GRCm39_Encode_Promoters.saf")
]

params.enhancer_saf_file = [
  "hs": file(params.GENOME_DIR + "/human/RegEmtDB_enhancer_hg38.saf"),
  "mm": file(params.GENOME_DIR + "/mouse/old_annots/GRCm39_Encode_Enhancers.saf")
]

params.alignment_ncore = 4
params.fragment_ncore = 6
params.ramsize = 5000000000
params.count_ncores = 3
params.fragment_ncores = 4


// subworkflows
include { 
  AlignPairedTag 
} from params.HOME_REPO + "/nf/workflows/v3/spikein_alignment.nf"
include { 
  PairedTagDNA as PrimaryDNA;
  PairedTagDNA as SpikeDNA
} from params.HOME_REPO + "/nf/workflows/v3/pairedtag_dna.nf"
include {
  PairedTagRNA as PrimaryRNA;
  PairedTagRNA as SpikeRNA
} from params.HOME_REPO + "/nf/workflows/v3/pairedtag_rna.nf"
include { 
  publishData as publish_rna_h5,
  publishData as publish_dna_h5
} from params.HOME_REPO + "/nf/modules/publish"

/* channel over rows of the digest */
pair_ch = Channel.fromPath(params.LIBRARY_DIGEST).splitCsv(header: true, sep: ",").map{ row -> tuple(row.sequence_id, file(row.fastq1), file(row.fastq2), row.lysis_id, row.library_type)}

workflow {
    aligned_bams = AlignPairedTag(pair_ch)
    prime_dna = PrimaryDNA(aligned_bams[0], params.SPECIES)
    spike_dna = SpikeDNA(aligned_bams[2], params.SPIKEIN_SPECIES)
    prime_rna = PrimaryRNA(aligned_bams[1], params.SPECIES)
    spike_rna = SpikeRNA(aligned_bams[3], params.SPECIES)
    publish_rna_h5(prime_rna[0])
    publish_dna_h5(prime_dna[0])
}
