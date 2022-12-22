"""
Merge tags from 2 bam files
"""
from argparse import ArgumentParser
import pysam


def get_args():
    parser = ArgumentParser('add_tags')
    parser.add_argument('barg1', help='String of the form <bamfile>:<source_tag>:<destination_tag>')
    parser.add_argument('barg2', help='String of the form <bamfile>:<source_tag>:<destination_tag>')
    parser.add_argument('out', help='the output bam file')
    parser.add_argument('--drop', help='comma-separated list of tags to drop', default=None)

    return parser.parse_args()


def main(args):
    bam1, src1, dest1 = args.barg1.split(':')
    bam2, src2, dest2 = args.barg2.split(':')
    drop = args.drop.split(',') if args.drop is not None else []
    b1h, b2h = pysam.AlignmentFile(bam1), pysam.AlignmentFile(bam2)
    hdr_dct = b1h.header.as_dict()
    out = pysam.AlignmentFile(args.out, mode='wb', header=hdr_dct)
    for r1, r2 in zip(b1h, b2h):
        assert r1.query_name == r2.query_name
        if r1.has_tag(src1):
            r1.set_tag(dest1, r1.get_tag(src1))
            r1.set_tag(src1, None)
        if r2.has_tag(src2):
            r1.set_tag(dest2, r2.get_tag(src2))
        for tg in drop:
            if r1.has_tag(tg):
                r1.set_tag(tg, None, replace=True)
        out.write(r1)

    out.close()


if __name__ == "__main__":
    main(get_args())

        
