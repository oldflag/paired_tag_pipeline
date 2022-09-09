"""
Streamlines the process of barcode parsing (R2 -> sample + well + well + UMI) and read-1 annotation
into a single process with no intermediate files.

The output annotated fastq files will be split by the sample ID

"""
from argparse import ArgumentParser
import parse_R2 as pr2
from itertools import islice
from multiprocessing import Pool
from utils import read_fasta, read_fastq, xopen, fastq_record


def get_args():
    parser = ArgumentParser('PairedTag_R2_Parser')
    parser.add_argument('R1_fastq', help='The Read 1 paired-tag fastQ, containing genomic sequence')
    parser.add_argument('R2_fastq', help='The R2 fastq file, containing barcodes')
    parser.add_argument('well_bc', help='The well barcode fasta file')
    parser.add_argument('sample_manifest', help='The sample manifest file, containing sample id, barcode, and antibody')
    parser.add_argument('linkers', help='The linker fasta file')
    parser.add_argument('outdir', help='The output directory')
    parser.add_argument('--threads', help='The number of threads to use', default=1, type=int)
    parser.add_argument('--umi_size', help='The size of the UMI in bp', default=10, type=int)
    parser.add_argument('--library_id', help='The library id to use (otherwise infer from fastq)', default=None, type=str)
    parser.add_argument('--sequence_id', help='The sequence id', default=None)
    parser.add_argument('--min_rlen', help='Minimum R1 read length for output', default=30, type=int)
    parser.add_argument('--umi_size', help='nt of UMI', default=10, type=int)

    return parser.parse_args()


def read_barcode_files(args):
    linker_seqs = list(read_fasta(args.linkers))
    combin_seqs = list(read_fasta(args.well_bc))
    sample_seqs = pr2.get_sample_seqs(args.sample_bc, args.R2_fastq, args.library_id)
    print(sample_seqs)

    pr2.check_args(args, linker_seqs, sample_seqs, combin_seqs)

    linker_size, sb_size, cb_size = (len(linker_seqs[0].seq),
                                     len(sample_seqs[0].seq),
                                     len(combin_seqs[0].seq))

    sample_seqmap = {sseq.seq: sseq.name for sseq in sample_seqs}
    combin_seqmap = {comb.seq: comb.name for comb in combin_seqs}
    return linker_seqs, combin_seqs, sample_seqs, linker_size, sb_size, cb_size, sample_seqmap, combin_seqmap


def iter_chunk_args(r1_fastq, r2_fastq, linker1, linker2, barcode_map,
                    sample_map, umi_size, chunk_size=10000):
    """
    For every `chunk_size` records in the fastq file, generate the arguments
    for a call to the R2 parse function.

    Inputs
    --------
    :r1_fastq: The path to the R1 fastq file
    :r2_fastq: The path to the R2 fastq file
    :linker1: The nucleotide sequence for linker1
    :linker2: The nucleotide sequence for linker2
    :barcode_map: A dictionary mapping sequence to labels (well barcodes)
    :sample_map: A dictionary mapping sequences to labels (sample barcodes)
    :umi_size: The size (in bp) of the UMI sequence
    :output_map: A dictionary of sample labels to out handles for writing annotated fastQ records

    Outputs
    ---------
    Iterator over argument tuples to the chunk parser
    """
    barcode_size = len(pr2.inext(barcode_map))
    sample_size = len(pr2.inext(sample_map))
    linker_size = len(linker1)
    r1_iter, r2_iter = read_fastq(r1_fastq), read_fastq(r2_fastq)
    r1_chunk = list(islice(r1_iter, chunk_size))
    r2_chunk = list(islice(r2_iter, chunk_size))
    while len(r1_chunk) > 0:
        args = (r1_chunk, r2_chunk, linker1, linker2, umi_size, barcode_size,
                sample_size, linker_size, barcode_map, sample_map)
        yield args
        r1_chunk = list(islice(r1_iter, chunk_size))
        r2_chunk = list(islice(r2_iter, chunk_size))


def annotate_reads(r1_recs, r2_recs, linker1, linker2, umi_size, barcode_size,
                   sample_size, linker_size, barcode_map, sample_map):
    barcodes = pr2.parse_barcodes_chunk(r2_recs, linker1, linker2, umi_size, barcode_size, sample_size,
                                        linker_size, barcode_map, sample_map)
    for read, (name, ident, barcodes) in zip(r1_recs, barcodes):
        rname = name + '|' + ident + '|' + barcodes
        if ident == '*':
            sm = '*'
        else:
            sm = ident.split(':')[-1]
        yield sm, fastq_record(rname, read.seq, read.qual)


def annotate_reads_(args):
    return annotate_reads(*args)


def main(args):
    # marshall the arguments
    linker_seqs, combin_seqs, sample_seqs, linker_size, sb_size, cb_size, sample_seqmap, combin_seqmap = read_barcode_files(args)
    # open a handle for each of the samples
    sequence_id = args.sequence_id if args.sequence_id is not None else args.R1_fastq.split('.')[0]
    annot_args = iter_chunk_args(args.R1_fastq, args.R2_fastq, linker_seqs[0], linker_seqs[1],
                                 combin_seqmap, sample_seqmap, args.umi_size)

    base = args.outdir + '/' + sequence_id
    handle_map = {sm: xopen(base + '_' + sm + '.fq.gz', 'wt')
                  for sm in sample_seqmap.values()}
    handle_map['*'] = xopen(base + '_unmatched.fq.gz', 'wt')

    if args.threads <= 1:
        chunk_iter = map(annotate_reads_, annot_args)
    else:
        pool = Pool(args.threads)
        chunk_iter = pool.imap(annotate_reads_, annot_args)

    for read_chunk in chunk_iter:
        for sample_id, record in read_chunk:
            handle_map[sample_id].write('@' + record.name + '\n')
            handle_map[sample_id].write(record.seq + '\n')
            handle_map[sample_id].write('+\n')
            handle_map[sample_id].write(record.qual + '\n')

    if args.threads > 1:
        pool.close()

    for sample_id, handle in handle_map.items():
        handle.close()


if __name__ == '__main__':
    main(get_args())
