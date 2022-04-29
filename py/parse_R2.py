"""
Parse R2 of a Paired-Tag experiment, given linker, well BC, and sample BC sequences.

The paired-tag experiment gives rise to two reads:

  R1:  Genomic sequence (DNA or RNA)
  R2:  Synthetic DNA used to determine cell and sample ID

The structure of R2 is as follows:

[UMI][BC1](0-3bp)[Linker1][BC2][Linker2](1bp)[SBC][RES](primer)

Where:
  UMI:  Unique Molecular Identifier (10bp)
  BC1:  Combinatorial barcode 1 (length determined by barcode design)
  Linker1: DNA linker sequence (30bp)
  BC2:  Combinatorial barcode 2 (length same as barcode 1)
  SBC:  Sample barcode (length determined by # of samples)
  RES:  Restriction enzyme binding site (8bp)
  Primer:  Any hanging primer sequences

The nature of the experiment is to tag cells with their sample
of origin (SBC), and then pool the cells into 96 or 384 well
plates, each of which contains its own barcode (BC1, BC2) which
is annealed to the DNA fragments. This step is done twice, generating
an ordered pair of barcodes (BC1, BC2) which (together with the
sample barcode) should disambiguate all cells.

This step 'Decodes' the cell coordinates and sample identity from the
sequence of Read2, by extracting the barcode sequences and cross
referencing them with the known barcode -> well and barcode -> sample
mappings (provided as fasta files).

"""
from argparse import ArgumentParser
from utils import read_fasta, read_fastq, xopen, fasta_record
from skbio.alignment import StripedSmithWaterman
from multiprocessing import Pool
from itertools import islice
from csv import DictReader

RE_LENGTH = 8

inext = lambda q: next(iter(q))

def get_args():
    parser = ArgumentParser('PairedTag_R2_Parser')
    parser.add_argument('R2_fastq', help='The R2 fastq file')
    parser.add_argument('well_bc', help='The well barcode fasta file')
    parser.add_argument('sample_bc', help='The sample barcode fasta file OR sample digest file (.csv)')
    parser.add_argument('linkers', help='The linker fasta file')
    parser.add_argument('output', help='The output .csv linking read name to disambiguated barcode')
    parser.add_argument('--threads', help='The number of threads to use', default=1, type=int)
    parser.add_argument('--umi_size', help='The size of the UMI in bp', default=10, type=int)

    return parser.parse_args()


def check_args(args, linker_seqs, sample_seqs, combin_seqs):
    # check that there are 2 linkers, and their size are equal
    if len(linker_seqs) > 2:
        raise ValueError('There should be only 2 linker sequences')

    if len(linker_seqs[0].seq) != len(linker_seqs[1].seq):
        raise ValueError('Size of linker sequences should match')

    if not all((len(sseq.seq) == len(sample_seqs[0].seq) 
                for sseq in sample_seqs)):
        raise ValueError('Length of sample barcode sequences must match')

    if not all((len(cseq.seq) == len(combin_seqs[0].seq)
                for cseq in combin_seqs)):
        raise ValueError('Length of combinatorial barcode sequences must match')


def iter_chunk_args(fastq_file, linker1, linker2, barcode_map,
                    sample_map, umi_size, chunk_size=5000):
    """
    For every `chunk_size` records in the fastq file, generate the arguments
    for a call to the R2 parse function.

    Inputs
    -------- 
    :fastq_file: The path to the R2 fastq file
    :linker1: The nucleotide sequence for linker1
    :linker2: The nucleotide sequence for linker2
    :barcode_map: A dictionary mapping sequence to labels (well barcodes)
    :sample_map: A dictionary mapping sequences to labels (sample barcodes)
    :umi_size: The size (in bp) of the UMI sequence

    Outputs
    ---------
    Iterator over argument tuples to the R2 chunk parser
    """
    barcode_size = len(inext(barcode_map))
    sample_size = len(inext(sample_map))
    linker_size = len(linker1)
    fastq_iter = read_fastq(fastq_file)
    chunk = list(islice(fastq_iter, chunk_size))
    while len(chunk) > 0:
        args = (chunk, linker1, linker2, umi_size, barcode_size,
                sample_size, linker_size, barcode_map, sample_map)
        yield args
        chunk = list(islice(fastq_iter, chunk_size))


def extract_barcodes(seq, linker1, linker2, lumi, lbc, lsn, lln):
    """
    Extract the barcodes from a Paired-Tag second read (R2).
    The format of this read is, loosely:

    [UMI][BC1](0-3bp)[Linker1][BC2][Linker2](1bp)[SBC][RES](primer)

    The relevant sequences can be extracted by the following process:

    i) Place Linker1 and Linker2 onto the sequence
    ii) Verify that Linker1 starts after BC1; otherwise
        left-pad the UMI with 'N'
    iii) Extract UMI from start of read
    iv) Extract BC1 from end of UMI
    v) Extract BC2 from start of Linker2
    vi) Extract SBC from (1bp after) Linker2 end
    vii) Extract RES from after SBC

    Input:
    -----------
    :seq:  The read sequence
    :linker1: a -StripedSmithWaterman- object for linker1
    :linker2: a -StripedSmithWaterman- object for linker2
    :lumi: the length of the UMI
    :lbc: the length of the (well) barcode
    :lsn: the length of the (sample) barcode
    :lln: the linker size


    Output:
    ----------
    :UMI: the UMI sequence
    :BC1: the BC1 sequence
    :BC2: the BC2 sequence
    :SBC: the sample barcode sequence
    :RES: the restriction enzyme sequence, pre-pended with
          the 1bp A/T base between Linker2 and SBC
    
    These values may be None under the following circumstances:
      + More than 3 'N' bases in read
      + Read length < umi_len + 2*bc_len + 2*linker_len + 1 + sample_len
      + Inability to properly place linker1 and linker2
      + linker 1 placed more than 5bp away from BC1

    """
    if len(seq) < (lumi + 2*lbc + 2*lln + lsn) or seq.count('N') > 3:
        return None, None, None, None, None

    l1 = linker1(seq)
    # force l2 to be aligned after l1
    seq2 = seq[(l1['target_begin'] + lln - 1):]
    l2 = linker2(seq2)
    l2 = {
        'optimal_alignment_score': l2['optimal_alignment_score'],
        'target_begin': l2['target_begin'] + lln + l1['target_begin'] - 1
    }

    if l2['target_begin'] - (l1['target_begin'] + lln) != lbc:
        # alignment failed - distance between linkers is not 1 barcode
        # if l2 alignment is good, just modify l1
        shift = lbc - (l2['target_begin'] - (l1['target_begin'] + lln))
        if (0 < shift < lbc) and l2['optimal_alignment_score'] > l1['optimal_alignment_score']:
            l1 = {'target_begin': l1['target_begin'] - shift}
        else:
            return None, None, None, None, None

    if l1['target_begin'] > lumi + lbc + 5:
        return None, None, None, None, None

    if l1['target_begin'] < lumi + lbc:
        # skipped base in UMI or barcode -- assume UMI
        bc1 = seq[(l1['target_begin']-lbc):l1['target_begin']]
        d_ = lumi + lbc - l1['target_begin']
        umi = 'N' * d_ + seq[:(lumi-d_)]
    else:
        umi, bc1 = seq[:lumi], seq[lumi:(lumi + lbc)]

    bc2 = seq[(l2['target_begin']-lbc):l2['target_begin']]
    eol2 = l2['target_begin'] + lln
    if eol2 + 1>= len(seq):
        # the read truncates on or after barcode 2
        tbase = 'N'
        sbc = 'N' 
        res = 'N' * RE_LENGTH
    else:
        tbase = seq[eol2]
        eol2 += 1
        sbc = seq[eol2:(eol2 + lsn)]
        res = seq[(eol2+lsn):(eol2+lsn+RE_LENGTH)]

    return umi, bc1, bc2, sbc, tbase + res


def parse_R2_barcodes_(args):
    return parse_barcodes_chunk(*args)

def parse_barcodes_chunk(reads, l1, l2, umi_bp, bc_bp, sn_bp, ln_bp,
                         bc_map, sn_map):
    """
    All R2 sequences in `reads` into UMI, well1, well2, sample, and
    restriction enzyme barcodes. See `parse_R2_barcodes`

    """
    swl1 = StripedSmithWaterman(l1)
    swl2 = StripedSmithWaterman(l2)
    return [parse_R2_barcodes(r, swl1, swl2, umi_bp, bc_bp, sn_bp,
                              ln_bp, bc_map, sn_map) for r in reads]


def parse_R2_barcodes(read, sw_l1, sw_l2, umi_bp, bc_bp, sn_bp, 
        ln_bp, bc_map, sn_map):
    """
    Parses R2 sequence into UMI, well1, well2, sample, and restriction
    enzyme barcodes, if possible.

    Input:
    ---------
    :read: a namedtuple containing the read name, sequence, and qualities
    :sw_l1: a -StripedSmithWaterman- object for linker1
    :sw_l2: a -StripedSmithWaterman- object for linker2
    :umi_bp: the size of the UMI (bp)
    :bc_bp: the size of the well barcodes (bp)
    :sn_bp: the size of the sample barcodes (bp)
    :ln_bp: the size of the linker sequences (bp)
    :bc_map: A dictionary mapping well barcodes to short IDs
    :sn_map: A dictionary mapping sample barcodes to IDs

    Output:
    -----------
    :arg1: the read name
    :arg2: the final barcode of the form {umi}:{bc1}:{bc2}:{sbc}:{dna/rna/unk};
           may be "*" for unknown or could not be mapped
    :arg3: the raw barcode; as above using raw sequences. May be "*" for 
           cases where the sequence could not be parsed.

    """
    umi, bc1, bc2, sbc, ers = \
            extract_barcodes(read.seq, sw_l1, sw_l2, umi_bp, bc_bp, sn_bp, ln_bp)

    if umi is None:
        return read.name, '*', '*'  # could not place linkers, or read filtered

    tr_umi = '*' if umi[0] == 'N' else umi
    tr_bc1 = bc_map[bc1] if bc1 in bc_map else rescue(bc1, bc_map)
    tr_bc2 = bc_map[bc2] if bc2 in bc_map else rescue(bc2, bc_map)
    tr_sbc = sn_map[sbc] if sbc in sn_map else rescue(sbc, sn_map)
    # the ers string consists of [typebp][RE seq]
    # if typebp is A, it is dna; T for rna; fallback is to use
    # the restriction enzyme sequence CCTGCAGG (dna) and GCGGCCGC (rna)
    tr_ers = 'rna' if 'TCGA' in ers else 'dna' if 'GGCC' in ers else 'unk'

    return (read.name, f'{tr_umi}:{tr_bc1}:{tr_bc2}:{tr_sbc}:{tr_ers}',
            f'{umi}:{bc1}:{bc2}:{sbc}:{ers}')


def rescue(seq, seqmap):
    """
    For a sequence not present in a sequence map `seqmap`, attempt to identify
    a corresponding sequence matching at all but 1 position, and return
    the value associated with that sequence

    """
    for seq2, val in seqmap.items():
        if sum((a != b for a, b in zip(seq, seq2))) <= 1:
            return val

    return '*'

def basename(fn):
    return fn.split('/')[-1]


def get_sample_seqs(fasta_or_digest, input_fastq):
    if fasta_or_digest[-4:] == '.csv':
        records = [r for r in DictReader(open(fasta_or_digest))]
        print(records)
        records = [r for r in records if basename(r['fastq2']) == basename(input_fastq)]
        print(records)
        return [fasta_record(r['assay_id'], r['barcode']) for r in records] 
    else:
        return list(read_fasta(fasta_or_digest))


def main(args):
    linker_seqs = list(read_fasta(args.linkers))
    combin_seqs = list(read_fasta(args.well_bc))
    sample_seqs = get_sample_seqs(args.sample_bc, args.R2_fastq)
    print(sample_seqs)

    check_args(args, linker_seqs, sample_seqs, combin_seqs)

    linker_size, sb_size, cb_size = (len(linker_seqs[0].seq),
                                     len(sample_seqs[0].seq),
                                     len(combin_seqs[0].seq))

    sample_seqmap = {sseq.seq: sseq.name for sseq in sample_seqs}
    combin_seqmap = {comb.seq: comb.name for comb in combin_seqs}

    arg_gen = iter_chunk_args(args.R2_fastq, linker_seqs[0].seq,
                linker_seqs[1].seq, combin_seqmap, sample_seqmap,
                args.umi_size, linker_size)

    if args.threads <= 1:
        demux_iter = map(parse_R2_barcodes_, arg_gen)
    else:
        pool = Pool(args.threads)
        demux_iter = pool.imap(parse_R2_barcodes_, arg_gen)

    with xopen(args.output, 'wt') as out:
        for elem in demux_iter:
            for readname, barcode, raw in elem:
                out.write(f'{readname},{barcode},{raw}\n')

    if args.threads > 1:
        pool.close()



if __name__ == '__main__':
    main(get_args())
