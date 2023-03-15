from argparse import ArgumentParser
import pysam
import re


def main(args):
    """
    Count the number and proportion of reads that contain the target
    tag; those that don't contain the target tag; those that contain
    the target tag and match the regex; and those that contain the
    target tag but don't match the regex
    """
    rex = re.compile(args.regex) if args.regex is not None else None
    counts = {'total': 0, 'Q30': 0, 'untagged': 0, 'tagged': 0, 'tagged_match': 0, 'tagged_nomatch': 0}
    bam = pysam.AlignmentFile(args.bam)
    for read in bam.fetch():
        counts['total'] += 1
        if read.mapping_quality >= args.min_mapq:
            counts['Q30'] += 1
            if read.has_tag(args.tag):
                counts['tagged'] += 1
                if rex is None or rex.match(read.get_tag(args.tag)) is not None:
                    counts['tagged_match'] += 1
                else:
                    counts['tagged_nomatch'] += 1
            else:
                counts['untagged'] += 1

    print('group\tcount\tproportion')
    for k in sorted(counts.keys()):
        print('\t'.join([k, str(counts[k]), '%.3f' % (counts[k]/counts['Q30'])]))


def get_args():
    parser = ArgumentParser()
    parser.add_argument('bam')
    parser.add_argument('tag')
    parser.add_argument('--regex', default=None)
    parser.add_argument('min_mapq', type=int, default=0)

    return parser.parse_args()


if __name__ == '__main__':
    main(get_args())
            
    

