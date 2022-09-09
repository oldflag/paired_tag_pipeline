"""
Filter and merge input bams on Q30 and no duplicate UMIs for MACS2-based peak calling.
UMIs will be filtered by retaining the left-most (first) read.

"""
from argparse import ArgumentParser
from add_tags import interval_key
import pysam


def read_key(rd, odr):
    return interval_key((rd.reference_name, rd.start, rd.end), odr)


def get_args():
    parser = ArgumentParser()
    parser.add_argument('bamlist')
    parser.add_argument('outbam')

    return parser.parse_args()


def extumi(rd):  # compute 'extended' umi (cell barcode + UMI)
    return rd.get_tag("CB") + ':' + rd.get_tag("MI", '.')


def umi_iter(bfile, wsize=300):
    hdl = pysam.AlignmentFile(bfile)
    stash, umis = list(), set()
    for read in hdl:
        if read.is_mapped and read.mapping_quality >= 30:
            umi = extumi(read)
            if umi not in umis:
                stash.append(read)
                umis.add(umi)
            i = 0
            for i, rd in enumerate(stash):
                if rd.reference_name != read.reference_name or read.start_position - rd.end_position >= wsize:
                    yield rd
                    umis.remove(extumi(rd))
                else:
                    break
            stash = stash[i:]
    for rd in stash:
        yield rd


class SortedIterator(object):
    def __init__(self, bam, iters):
        self.codr_ = {k: i for i, k in enumerate(bam.references)}
        self.iters = list()
        self.heads = list()
        self.vals = list()
        for itr in iters:
            self.add_iterator(itr)

    def add_iterator(self, iterator):
        x = next(iterator)
        self.iters.append(iterator)
        self.heads.append(x)
        self.vals.append(read_key(x, self.codr_))

    def __next__(self):
        ix, vl = min(enumerate(self.vals), key=lambda p: p[1])
        yield self.heads[ix]
        try:
            self.heads[ix] = next(self.iters[ix])
            self.vals[ix] = read_key(self.heads[ix])
        except StopIteration:
            del self.heads[ix]
            del self.vals[ix]
            del self.iters[ix]

    def __iter__(self):
        return self


def main(args):
    bams = [umi_iter(line.strip()) for line in open(args.bamlist)]
    firstbam = pysam.AlignmentFile(next(open(args.bamlist)).strip())
    sort_iter = SortedIterator(firstbam, bams)

    out_hdl = pysam.AlignmentFile(args.outbam, mode='wb', template=firstbam)
    firstbam.close()
    for read in iter(sort_iter):
        out_hdl.write(read)

    out_hdl.close()
