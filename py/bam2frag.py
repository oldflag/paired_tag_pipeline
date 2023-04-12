# Convert a bam file into a 10x-style fragment (bed) file
# this is to be run on the DNA libraries where UMI co-occur at
# the same location
from argparse import ArgumentParser
from multiprocessing import Pool
import pysam
import numpy as np
import sys
from collections import OrderedDict, Counter, defaultdict
import os
import gzip

FRAGMENT_SIZE = 200

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
    umi_positions = defaultdict(list)
    lastpos, nproc = start, 1
    cell_counts, cell_reads = Counter(), Counter()
    for read in fetch(bampath, (contig, start, end)):
       if read.is_mapped and read.has_tag('MI') and (read.reference_end - read.reference_start) >= 20:
           cell = read.get_tag('CB')
           mi = read.get_tag('MI')
           umi=(cell, mi)
           umi_positions[umi].append(read.reference_start)
           if read.reference_start > lastpos:
               nproc += 1
               lastpos = read.reference_start
               if nproc >= FRAGMENT_SIZE:  # we have moved a fragment, clean the dictionary
                   nproc = 1
                   out_umi = list()
                   for ukey, poslist in list(umi_positions.items()):
                       if read.reference_start - poslist[-1] > FRAGMENT_SIZE:
                           # it's been more than 1 fragment since we last saw this; count it
                           frag_start = int(sum(poslist)/len(poslist))
                           out_umi.append((contig, frag_start, frag_start + FRAGMENT_SIZE, ukey[0], len(poslist)))
                           del umi_positions[ukey]
                       elif read.reference_start - poslist[0] > FRAGMENT_SIZE:
                           # weird edge case - the UMI has been seen twice with about a fragment length
                           # difference or more; cut the position list and count it
                           seg1, seg2 = list(), list()
                           for p in poslist:
                               if (p - poslist[0]) >= (poslist[-1] - p):
                                   seg1.append(p)
                               else:
                                   seg2.append(p)
                           frag_start = int(sum(seg1)/len(seg1))
                           out_umi.append((contig, frag_start, frag_start + FRAGMENT_SIZE, umi[0], len(seg1)))
                           umi_positions[ukey] = seg2 
                   out_umi.sort(key=lambda k: k[1])
                   for orec in out_umi:
                       cell_counts[orec[-2]] += 1
                       cell_reads[orec[-2]] += orec[-1]
                       ohdl.write('\t'.join(map(str, orec)) + '\n')


    for umi, poslist in list(umi_positions.items()):
        frag_start = int(sum(poslist)/len(poslist))
        record = (contig, frag_start, FRAGMENT_SIZE + frag_start, umi[0], len(poslist))
        cell_counts[record[-2]] += 1
        cell_reads[record[-2]] += record[-1]
        ohdl.write('\t'.join(map(str, record)) + '\n')


    ohdl.close()
    return outgz, cell_counts, cell_reads


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

    tot_counts, tot_umi = dict(), dict()
    for _, cell_umi, cell_reads in tx_chunks:
        for cell, cnt in cell_umi.items():
            tot_counts[cell] = tot_counts.get(cell, 0) + cell_reads[cell]
            tot_umi[cell] = tot_umi.get(cell, 0) + cnt

    nr = sum(tot_counts.values())
    nu = sum(tot_umi.values())
    print('%s: Counted %d aligned, barcoded reads supporting %d umi (efficiency %.1f%%)' % (args.bam, nr, nu, 100*nu/nr))

    valid_bcs = {k for k, v in tot_umi.items() if v >= args.min_count}
    if known_bcs is not None or args.min_count > 1:
        valid_bcs = valid_bcs & (known_bcs or valid_bcs)
        print("Reduced %d observed barcodes to %d filtered" % (len(tot_counts), len(valid_bcs)))


    rmlist=list()

    hdl = gzip.open(args.tsvgz, 'wt') 
    for f_chunk, _, _ in tx_chunks:
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
