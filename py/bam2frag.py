# Convert a bam file into a 10x-style fragment (bed) file
# this is to be run on the DNA libraries where UMI co-occur at
# the same location
from argparse import ArgumentParser
from multiprocessing import Pool
import pysam
import numpy as np
import sys
from collections import OrderedDict
import os
import gzip

FRAGMENT_SIZE = 100

def get_args():
    parser = ArgumentParser()
    parser.add_argument('bam', help='The input bam file')
    parser.add_argument('tsvgz', help='The output tsv.gz file')
    parser.add_argument('--ncores', help='The number of cores to use', default=1, type=int)

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
    cell_umi_locus = OrderedDict()
    for read in fetch(bampath, (contig, start, end)):
       if read.has_tag('MI'):
           if read.is_forward:
               read_loc = (read.reference_name, read.reference_start, read.reference_start + FRAGMENT_SIZE)
           else:
               read_loc = (read.reference_name, max(1, read.reference_end - FRAGMENT_SIZE), read.reference_end)
           if loc[0] != read_loc[0] or loc[1] != read_loc[1]:
               for (cell, loc, umi), count  in cell_umi_locus.items():
                   ohdl.write('%s\t%d\t%d\t%s\t%d\n' % (loc[0], loc[1], loc[2], cell, count))
               loc, cell_umi_locus = read_loc, OrderedDict()
           umi = read.get_tag('MI')
           cell = read.get_tag('CB')
           cell_umi_locus[(cell, read_loc, umi)]  = cell_umi_locus.get((cell, read_loc, umi),0) + 1

    for (cell, loc, umi), count in cell_umi_locus.items():
        ohdl.write('%s\t%d\t%d\t%s\t%d\n' % (loc[0], loc[1], loc[2], cell, count))

    ohdl.close()
    return outgz


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
    chunk_iter = chunk_reads(args.bam, args.tsvgz, 25000000) # 250M bp/chunk
    if args.ncores > 1:
        pool = Pool(args.ncores)
        tx_chunks = pool.map(chunk_to_frags, chunk_iter)
    else:
        tx_chunks = map(chunk_to_frags, chunk_iter)

    rmlist=list()
    hdl = gzip.open(args.tsvgz, 'wt') 
    for f_chunk in tx_chunks:
        print(f_chunk)
        in_ = gzip.open(f_chunk, 'rt')
        hdl.writelines(in_)
        in_.close()
        rmlist.append(f_chunk)
    hdl.close()

    os.system('rm %s' % ' '.join(rmlist))


if __name__ == '__main__':
    main(get_args())
