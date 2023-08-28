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
  merge_counts as merge_bin_read;
  merge_counts as merge_peak_read
} from params.HOME_REPO + "/nf/modules/count_v3"

include {
  publishData as publishtracks
  catAndPublish as publishdup_dna
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
       suffix // any suffix (such as spike-in)
  main:
       params.macs_genome_type = genome
       // group by antibody
       antibody_grp = dna_ch.map{ it -> tuple(it[4], it[1])}.filter{ ! it[0].contains("UNK")}.groupTuple().map{ it -> tuple(it[0], it[1], params.RUN_NAME)}

       // call peaks
       peaks = MACS2_multi(antibody_grp, genome, params.genome_reference_dir[genome],params.genome_reference_name[genome], params.py_dir, params.sh_dir )  // input: (antibody_name, bam_file_list, experiment_name)

       // merge peaks prior to tagging
       mg_peaks = merge_saf(peaks[0].map{ it -> it[2]}.collect(), 'all_antibodies_' + genome)

       // tag the original dna bams
       tagged_ch = add_tags(dna_ch.filter{ ! it[1].simpleName.contains("_unlinked") },
                            tuple(params.genome_bin_file[genome], "BN", "SAF"),
                            tuple(params.enhancer_saf_file[genome], "EE", "SAF"),
                            tuple(params.promoter_saf_file[genome], "EP", "SAF"),
                            tuple(params.genome_saf_file[genome], "GN", "SAF"),
                            mg_peaks.map{it -> tuple(it, "PK", "SAF")},
                            tuple(file("input.6"), "NA", "NAN"),
                            params.py_dir,
                            params.sh_dir)

       publishdup_dna(tagged_ch[1].map{ it[1]}.collect(), params.RUN_NAME + suffix + '.duplication_metrics.dna.csv')

       // strategy:  (antibody, assay, alignment, bam).join(antibody, saf).map(assay"_"antibody, alignment, bam, saf)
       qc_input = inner_join(tagged_ch[0].map{ it -> tuple(it[4], it[3], it[0], it[1]) },
                             peaks[0].map{ it -> tuple(it[3], it[2]) }).map{
                                 it -> tuple(it[2] + '_' + it[0] + '_' + it[1], it[3], it[4])}
 
       qc_output = chip_qc(qc_input, params.py_dir)
 
       // merge the files
       merge_chip_qc(qc_output.map{it -> it[2]}.collect(), 
                     qc_output.map{it -> it[1]}.collect(), 
                     params.RUN_NAME + "_" + genome,
                     params.py_dir)
 
       // these have both umi and read counts
       bin_counts = dna_count(tagged_ch[0], "BN")  
       peak_counts = peak_count(tagged_ch[0], "PK")
 
       bin_h5ad =  merge_bin(prep_h5_inputs(bin_counts[0], 3, "DNA_Q30_UMICount_per_bin_" + genome),
                             params.SAMPLE_DIGEST,params.py_dir)
       peak_h5ad = merge_peak(prep_h5_inputs(peak_counts[0], 3, "DNA_Q30_UMICount_per_peak_" + genome),
                              params.SAMPLE_DIGEST, params.py_dir)

       bin_read_h5 = merge_bin_read(prep_h5_inputs(bin_counts[0], 4, "DNA_Q30_READCount_per_bin_" + genome),
                                    params.SAMPLE_DIGEST, params.py_dir)
       peak_read_h5 = merge_peak_read(prep_h5_inputs(peak_counts[0], 4, "DNA_Q30_READCount_per_peak_" + genome),
                                      params.SAMPLE_DIGEST, params.py_dir)
    
  emit:
      peak_umi_h5ad = peak_h5ad[0]
      bin_umi_h5ad = bin_h5ad[0]
      bin_read_h5ad = bin_read_h5[0]
      peak_umi_h5ad = peak_read_h5[0]
      
}
