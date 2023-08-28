nextflow.enable.dsl=2

include { 
  add_tags  
} from params.HOME_REPO + "/nf/modules/pairedtag_reads_v3"

include { 
  umitools_count as rna_count;
  merge_counts as rna_merge_counts; 
  merge_counts as rna_merge_readcounts;
} from params.HOME_REPO + "/nf/modules/count"

include {
  catAndPublish as publishdup_rna
} from params.HOME_REPO + "/nf/modules/publish"


workflow PairedTagRNA {
  take:
    rna_ch   // channel over rna bams
    genome   // genome abbrev ("hs" // params.SPECIES)
    suffix   // any suffix for the call, such as "spikein"

  main:

    rna_tagged = add_tags(rna_ch.filter{ ! it[1].simpleName.contains("_unlinked") },
                         tuple(params.genome_bin_file[genome], "BN", "SAF"),
                         tuple(params.genome_saf_file[genome], "GN", "SAF"),
                         tuple(file("input.3"), "NA", "NAN"),
                         tuple(file("input.4"), "NA", "NAN"),
                         tuple(file("input.5"), "NA", "NAN"),
                         tuple(file("input.6"), "NA", "NAN"),
                         params.py_dir,
                         params.sh_dir)

    publishdup_rna(rna_tagged[1].map{ it[1]}.collect(), params.RUN_NAME + suffix + '.duplication_metrics.rna.csv')

  
    // read and umi count with umi_tools based on a given tag //
    // RNA read and umi count per cell per gene
    rna_counts = rna_count(rna_tagged[0], "GN")
    r_umi_input = rna_counts[0].map{it -> 
        tuple(1, it[0], it[1], it[2], it[3])
     }.groupTuple().map{ 
        it -> tuple(it[3].collect(), "RNA_Q30_UMICount_per_gene_" + genome)
    }
    rna_umi_merged_h5ad = rna_merge_counts(r_umi_input, params.SAMPLE_DIGEST, params.py_dir)
    r_read_input = rna_counts[0].map{it -> 
        tuple(1, it[0], it[1], it[2], it[3])
     }.groupTuple().map{ 
        it -> tuple(it[4].collect(), "RNA_Q30_READCount_per_gene_" + genome)
    }
    rna_read_merged_h5ad = rna_merge_readcounts(r_read_input, params.SAMPLE_DIGEST, params.py_dir)
  
  emit:
    rna_umi_h5ad = rna_umi_merged_h5ad[0]
    rna_read_h5ad = rna_read_merged_h5ad[0]

}
  
