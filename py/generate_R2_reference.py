"""
Builds a fasta reference file from user-provided barcodes. The inputs are:
    1) A fasta file containing valid barcodes for 96/384 well plates
    2) A fasta file containing valid barcodes for individual samples
    3) A fasta file containing linker sequences

The reference will be generated for all combinatorial sequences of the form:
    ([L1])[WELL_BC][L2][WELL_BC][L3][SAMPLE_BC]

Note that Linker 1 is optional -- if 2 linkers are provided there will
be no prefix linker; if 3 are provided there will be a linker sequence.

Note that the read itself is prefixed by a UMI; and has a 2bp suffix
from the restriction enzyme site. These are to be stripped prior to
alignment.

"""
from argparse import ArgumentParser
from collections import namedtuple
from utils import read_fasta
from itertools import product as iproduct

def get_args():
    parser = ArgumentParser('Build a reference for R2 from barcodes & linkers\npython generate_R2_reference.py ')
    parser.add_argument('well_barcode_fasta', help='FASTA file containing valid well barcode sequences')
    parser.add_argument('sample_barcode_fasta', help='FASTA file containing valid sample barcode sequences')
    parser.add_argument('--linkers', help='FASTA file containing linker sequences (optional)', default=None)

    return parser.parse_args()


def main(args):
    wells = list(read_fasta(args.well_barcode_fasta))
    samples = list(read_fasta(args.sample_barcode_fasta))
    if args.linkers:
        linkers = list(read_fasta(args.linkers))
    else:
        linkers = []

    if len(linkers) == 3:
        pfx, linkers = linkers[0], linkers[1:]
    else:
        pfx = namedtuple('whatever', ('seq'))('')

    for w1, w2, s in iproduct(wells, wells, samples):
        bcname = '>%s:%s:%s' % (w1.name, w2.name, s.name)
        if linkers:
            bcseq = f'{pfx.seq}{w1.seq}{linkers[0].seq}{w2.seq}{linkers[1].seq}{s.seq}'
        else:
            bcseq = f'{w1.seq}{w2.seq}{s.seq}'
        print(bcname)
        print(bcseq)


if __name__ == '__main__':
    main(get_args())
