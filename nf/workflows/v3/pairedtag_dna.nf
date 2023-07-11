nextflow.enable.dsl=2

// modules
include {
  MACS2_multi;
  merge_saf;
  chip_qc;
  merge_chip_qc
} from params.HOME_REPO + "/nf/modules/peaks"

include {
  add_tags
} from params.HOME_REPO + "/nf/modules/pairedtag_reads_v3"

include {
  umitools_count as dna_count;
  umitools_count as peak_count;
  merge_counts as merge_bin;
  merge_counts as merge_peak;
} from params.HOME_REPO + "/nf/modules/count_v3"

include {
  publishData as publishtracks
} from params.HOME_REPO + "/nf/modules/publish"


def inner_join(ch_a, ch_b) {
    return ch_b.cross(ch_a).map { [it[0][0], *it[1][1..-1], *it[0][1..-1]] }
}


def prep_h5_inputs(count_ch, idx, name) {
    return count_ch.map{
        it -> tuple(1, it[0], it[1], it[2], it[3])
    }.groupTuple().map{
        it -> tuple(it[idx].collect(), name)
    }
}

workflow PairedTagDNA {
  take:
       dna_ch  // dna over dna bam files
       genome // the genome abbrev (params.SPECIES)
  main:
       params.macs_genome_type = genome
       // group by antibody
       antibody_grp = dna_ch.map{ it -> tuple(it[4], it[1])}.groupTuple().map{ it -> tuple(it[0], it[1], params.RUN_NAME)}

       // call peaks
       peaks = MACS2_multi(antibody_grp, genome)  // input: (antibody_name, bam_file_list, experiment_name)

       // merge peaks prior to tagging
       mg_peaks = merge_saf(peaks[0].map{ it -> it[2]}.collect(), 'all_antibodies_' + genome)

       // tag the original dna bams
       tagged_ch = add_tags(dna_ch.filter{ ! it[1].simpleName.contains("_unlinked") },
                            tuple(params.genome_bin_file[genome], "BN", "SAF"),
                            tuple(params.enhancer_saf_file[genome], "EE", "SAF"),
                            tuple(params.promoter_saf_file[genome], "EP", "SAF"),
                            tuple(params.genome_saf_file[genome], "GN", "SAF"),
                            mg_peaks.map{it -> tuple(it, "PK", "SAF")},
                            tuple(file("input.6"), "NA", "NAN"))

       // strategy:  (antibody, assay, alignment, bam).join(antibody, saf).map(assay"_"antibody, alignment, bam, saf)
       qc_input = inner_join(tagged_ch.map{ it -> tuple(it[4], it[3], it[0], it[1]) },
                             peaks[0].map{ it -> tuple(it[3], it[2]) }).map{
                                 it -> tuple(it[2] + '_' + it[0] + '_' + it[1], it[3], it[4])}
 
       qc_output = chip_qc(qc_input)
 
       // merge the files
       merge_chip_qc(qc_output.map{it -> it[2]}.collect(), 
                     qc_output.map{it -> it[1]}.collect(), 
                     params.RUN_NAME + "_" + genome)
 
       bin_counts = dna_count(tagged_ch, "BN")
       peak_counts = peak_count(tagged_ch, "PK")
 
       bin_h5ad =  merge_bin(prep_h5_inputs(bin_counts[0], 3, "DNA_Q30_UMICount_per_bin_" + genome),
                             params.SAMPLE_DIGEST)
       peak_h5ad = merge_peak(prep_h5_inputs(peak_counts[0], 3, "DNA_Q30_UMICount_per_peak_" + genome),
                              params.SAMPLE_DIGEST)
    
  emit:
      peaks = peak_h5ad[0]
      bins = bin_h5ad[0]
}
