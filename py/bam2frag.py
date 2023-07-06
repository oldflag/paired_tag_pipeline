# Convert a bam file into a 10x-style fragment (bed) file
# this is to be run on the DNA libraries where UMI co-occur at
# the same location
from argparse import ArgumentParser
from multiprocessing import Pool
import pysam
import numpy as np
import sys
from collections import OrderedDict, Counter
import os
import gzip

FRAGMENT_SIZE = 100

def get_args():
    parser = ArgumentParser()
    parser.add_argument('bam', help='The input bam file')
    parser.add_argument('tsvgz', help='The output tsv.gz file')
    parser.add_argument('--ncores', help='The number of cores to use', default=1, type=int)
    parser.add_argument('--barcodes', help='Barcodes to retain', default=None)
    parser.add_argument('--min_count', help='The minimum UMI for filtered barcodes', default=0, type=int)

    return parser.parse_args()



def fetch(bamf, loc):
    hdl = pysam.AlignmentFile(bamf)
    itr = hdl.fetch(loc[0], loc[1], loc[2])
    # do not include reads starting before the start
    try:
         read = next(itr)
    except StopIteration:
         read = None
    while read and read.reference_start < loc[1]:
        try:
            read = next(itr)
        except StopIteration:
            read = None
            break
    if read:
        yield read
    for read in itr:
        yield read
    

def chunk_to_frags(bam_with_chunk):
    bampath, contig, start, end, outgz = bam_with_chunk
    ohdl = gzip.open(outgz, mode='wt')
    loc = ('None', 0, 0)
    cell_umi_locus, cell_counts = OrderedDict(), Counter()
    for read in fetch(bampath, (contig, start, end)):
       if read.has_tag('MI'):
           #if read.is_forward:
           #    read_loc = (read.reference_name, read.reference_start, read.reference_start + FRAGMENT_SIZE)
           #else:
           #    read_loc = (read.reference_name, max(1, read.reference_end - FRAGMENT_SIZE), read.reference_end)
           read_loc = (read.reference_name, read.reference_start, read.reference_end)
           if read_loc[2] - read_loc[1] < 20:
               continue
           if loc[0] != read_loc[0] or loc[1] != read_loc[1]:
               for (cell, loc, umi), count  in cell_umi_locus.items():
                   ohdl.write('%s\t%d\t%d\t%s\t%d\n' % (loc[0], loc[1], loc[2], cell, count))
               loc, cell_umi_locus = read_loc, OrderedDict()
           umi = read.get_tag('MI')
           cell = read.get_tag('CB')
           cell_umi_locus[(cell, read_loc, umi)]  = cell_umi_locus.get((cell, read_loc, umi),0) + 1
           cell_counts[cell] += 1

    for (cell, loc, umi), count in cell_umi_locus.items():
        ohdl.write('%s\t%d\t%d\t%s\t%d\n' % (loc[0], loc[1], loc[2], cell, count))

    ohdl.close()
    return outgz, cell_counts


def chunk_reads(bam, base_file, chunk_size):
    hdl = pysam.AlignmentFile(bam)
    contigs = hdl.references
    lens = hdl.lengths
    hdl.close()

    for contig, length in zip(contigs, lens):
        breaks = [x for x in np.arange(start=0, stop=length, step=chunk_size)] + [length]
        starts = breaks[:-1]
        ends = breaks[1:]
        for start, end in zip(starts, ends):
            if start >= end:
                print('bad: %s:%d-%d' % (contig, start, end))
                continue
            of = base_file[:-len('.tsv.gz')] + '_{}_{}_{}.tsv.gz'.format(contig, start, end)
            yield (bam, contig, start, end, of)


def main(args):
    chunk_iter = chunk_reads(args.bam, args.tsvgz, 250000000) # 250M bp/chunk
    if args.ncores > 1:
        pool = Pool(args.ncores)
        tx_chunks = pool.map(chunk_to_frags, chunk_iter)
    else:
        tx_chunks = map(chunk_to_frags, chunk_iter)

    if args.barcodes is not None:
        known_bcs = {x.strip() for x in open(args.barcodes)}
    else:
        known_bcs = None

    tot_counts = dict()
    for _, ccounts in tx_chunks:
        for cell, cnt in ccounts.items():
            tot_counts[cell] = tot_counts.get(cell, 0) + cnt

    valid_bcs = {k for k, v in tot_counts.items() if v >= args.min_count}
    if known_bcs is not None:
        valid_bcs = valid_bcs & known_bcs

    print("Reduced %d observed barcodes to %d filtered" % (len(tot_counts), len(valid_bcs)))

    rmlist=list()

    hdl = gzip.open(args.tsvgz, 'wt') 
    for f_chunk, _ in tx_chunks:
        in_ = gzip.open(f_chunk, 'rt')
        for line_ in in_:
            fields = line_.strip().split('\t')
            if fields[3] in valid_bcs:
                hdl.write(line_)
        in_.close()
        rmlist.append(f_chunk)
    hdl.close()

    os.system('rm %s' % ' '.join(rmlist))


if __name__ == '__main__':
    main(get_args())
