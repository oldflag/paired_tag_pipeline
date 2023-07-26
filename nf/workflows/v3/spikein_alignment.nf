nextflow.enable.dsl=2

// modules
include { trim_fq_single as trim_rna; trim_fq_single as trim_dna} from params.HOME_REPO + "/nf/modules/trim"
include { process_pairedtag; 
          barcode_qc as barcode_qc} from params.HOME_REPO + '/nf/modules/pairedtag_reads_v3'
include { star_aligner; bwa_aligner; bam_to_frag} from params.HOME_REPO + "/nf/modules/alignment_spikein"
include { alignment_qc as alignment_qc_spike; 
          alignment_qc as alignment_qc_primary; 
          merge_alignment_qc as merge_alignment_qc_spike;
          merge_alignment_qc as merge_alignment_qc_primary } from params.HOME_REPO + "/nf/modules/alignment"
include { rnaseqc_call as rnaseqc_primary; 
          rnaseqc_call as rnaseqc_spike;
          merge_rnaseqc } from params.HOME_REPO + "/nf/modules/rnaseqc_spikein"
include { catAndPublishDF as publishspikeqc;
          publishData as publishfragments;
          publishData as publishlogo;
          catAndPublish as publishlogs;
          catAndPublish as publishfraghist;
          publishData as publishrnaqc;
          publishData as publishbarcodeqc;
          publishData as publishalignmentqc;
          publishData as publishchipqc;
          publishData as publishtracks } from params.HOME_REPO + "/nf/modules/publish" 

/* channel over rows of the digest */


workflow AlignPairedTag {
  take:
      pair_ch  // channel over read pairs: (seqid, fq1, fq2, lysid, type)
  main:
      type_ch = pair_ch.map{ it -> tuple(it[0], it[4])}
      split_fqs = process_pairedtag(pair_ch.map{ it -> tuple(it[0], it[1], it[2], it[3])}, params.py_dir, params.combin_barcodes, params.sample_barcodes, params.linker_file)
      publishlogo(split_fqs[3])
      i=0
      j=0
      k=0
      fastq_ids = split_fqs[0].map{ it -> [i++, it] }
      fastq_r1files = split_fqs[1].map{ it -> [j++, it]}
      fastq_r2files = split_fqs[2].map{ it -> [k++, it]}
       // join them
      fqjoin = fastq_ids.join(fastq_r1files).map{ it -> tuple(it[0], it[1], it[2])}
      fqjoin = fqjoin.join(fastq_r2files).map{ it -> tuple(it[1], it[2], it[3]) }
      fqjoin = fqjoin.join(type_ch).map{ it -> tuple(it[0], it[3], it[1], it[2]) }
      fqjoin = fqjoin.transpose()
      // fqjoin.subscribe{ println it }

      // this is now a tuple of (seq_id, seq_type, fastq1, fastq2)
     
      barcode_pdfs = barcode_qc(fqjoin.map{ it -> tuple(it[0], it[2]) }.groupTuple(), params.py_dir, params.plate_layout)
      publishbarcodeqc(barcode_pdfs)
      
      dna_fq = fqjoin.filter{ it[1] =~ /dna/ }.map{it -> tuple(it[0], it[2], it[3], it[1])}
      dna_raw_bams = bwa_aligner(dna_fq, params.bwa_index[params.SPECIES], params.genome_reference[params.SPECIES], params.bwa_index[params.SPIKEIN_SPECIES],params.genome_reference[params.SPIKEIN_SPECIES], params.py_dir)
      
      rna_fq = fqjoin.filter{ it[1] =~ /rna/ }.map{it -> tuple(it[0], it[2], it[3], it[1])}
      //rna_fq.subscribe{ println(it) }
      rna_raw_bams = star_aligner(rna_fq, params.star_index[params.SPECIES], params.star_index[params.SPIKEIN_SPECIES], params.py_dir)

      //rna_qc
      rnaqc_primary = rnaseqc_primary(rna_raw_bams[0].map{it -> tuple(it[0], it[1], it[3], it[4])},
                                      params.genome_gtf_collapsed_file[params.SPECIES])
      rnaqc_spike = rnaseqc_spike(rna_raw_bams[1].map{it -> tuple(it[0], it[1], it[3], it[4])},
                                  params.genome_gtf_collapsed_file[params.SPIKEIN_SPECIES])
      rnaqc_mg = merge_rnaseqc(rnaqc_primary.mix(rnaqc_spike).map{it -> it[4]}.collect(), params.RUN_NAME, params.py_dir)
    
      publishrnaqc(rnaqc_mg[1])
    
      // in v2 pipeline, a single .bam is produced
      // in v3 pipeline, a (filtered) primary bam is produced in slot 0
      //                 a (filtered) spike-in bam i sproduced in slot 1
      //                 a spike-in info file is produced in slot 2
    
      publishspikeqc(dna_raw_bams[2].mix(rna_raw_bams[2]).collect(), params.RUN_NAME + "_spikein_assignment_qc.csv")
    
      primary_bams_all = dna_raw_bams[0].mix(rna_raw_bams[0])
      spikein_bams_all = dna_raw_bams[1].mix(rna_raw_bams[1])
    
      bamqc = alignment_qc_primary(primary_bams_all)
      bamqc_spike = alignment_qc_spike(spikein_bams_all)
      alignment_qcfile = merge_alignment_qc_primary(bamqc.map{ it -> it[1]}.collect(), params.RUN_NAME)
      alignment_qcfile_spike = merge_alignment_qc_spike(bamqc_spike.map{ it -> it[1]}.collect(), params.RUN_NAME + '_spikein')
      publishalignmentqc(alignment_qcfile.mix(alignment_qcfile_spike))
    
      // fragments
      fragfiles = bam_to_frag(dna_raw_bams[0].mix(dna_raw_bams[1]).map{ it -> it[1] }, params.py_dir)  // [0]: fragments [1]: index [2]: logs
      publishfragments(fragfiles[0])
      publishfraghist(fragfiles[4].collect(), params.RUN_NAME + '_fragment_sizes.txt')
      publishlogs(fragfiles[2].collect(), "gather_convert_fragments.log")
    
    
  emit:
      primary_dna = dna_raw_bams[0]
      primary_rna = rna_raw_bams[0]
      spikein_dna = dna_raw_bams[1]
      spikein_rna = rna_raw_bams[1]
}
