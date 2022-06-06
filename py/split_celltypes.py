# Convert a bam file into a 10x-style fragment (bed) file
# this is to be run on the DNA libraries where UMI co-occur at
# the same location
from argparse import ArgumentParser
from multiprocessing import Pool
import pysam
import sys
from collections import OrderedDict
import os
import gzip


def get_args():
    parser = ArgumentParser()
    parser.add_argument('bam', help='The input bam file')
    parser.add_argument('celltype_map', help='A text file of the form <cell_barcode>,<cell_type>')
    parser.add_argument('--ncores', help='The number of cores to use', default=1, type=int)

    return parser.parse_args()


def subset_type(cell_type, type_map, bam_path):
    in_bam = pysam.AlignmentFile(bam_path, 'r')
    cell_map = ok_cells = {k for k, v in type_map.items() if v == cell_type}
    out_path = bam_path[:-4] + '_' + cell_type + '.bam'
    in_hdr = in_bam.header
    cbc1, cbc2 = set(), set()
    with pysam.AlignmentFile(out_path, 'wb', header=in_hdr) as outr:
        for i, rd in enumerate(in_bam.fetch()):
            bc = rd.get_tag('CB')
            if bc not in cell_map:
                continue
            cbc = str(rd.reference_start) + bc
            if cbc in cbc1 or cbc in cbc2:
                continue
            outr.write(rd)
            if (i+1) % 1000000 == 0:
                cbc2, cbc1 = cbc1, set()  # clear the back half of the cache
    in_bam.close()
    return True 
             

def subset_type_(args):
    return subset_type(*args)


def main(args):
    types, typelist = dict(), set()
    for entry in open(args.celltype_map, 'rt'):
        cell, type_ = entry.strip().split(',')
        types[cell] = type_
        typelist.add(type_)

    if args.ncores > 1:
        pool = Pool(args.ncores)
        success = pool.map(subset_type_, ((type_, types, args.bam) for type_ in typelist))
    else:
        success = map(subset_type_, ((type_, types, args.bam) for type_ in typelist))

    pool.close()


if __name__ == '__main__':
    main(get_args())
