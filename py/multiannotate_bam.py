from argparse import ArgumentParser
from collections import OrderedDict
import pysam


def interval_key(loc1, ctg_odr):
    return 10**9 * ctg_odr[loc1[0]] + loc1[1] + 10**(-6) * (loc1[2] - loc1[1])


def overlaps(loc1, loc2, ctg_odr=None):
    return (loc1[0] == loc2[0]) and (loc2[2] >= loc1[1]) and (loc2[1] <= loc1[2])


def is_before(loc1, loc2, ctg_odr):
    return ctg_odr[loc1[0]] < ctg_odr[loc2[0]] or (loc1[0] == loc2[0] and loc1[2] < loc2[1])


def is_after(loc1, loc2, ctg_odr):
    return ctg_odr[loc1[0]] > ctg_odr[loc2[0]] or (loc1[0] == loc2[0] and loc1[1] > loc2[2])


def load_track(interval_file, format, ref_order):
    if format == 'BED':
        ix_ = [0, 1, 2, 3]
    elif format == 'SAF':
        ix_ = [1, 2, 3, 0]
    else:
        raise ValueError(f'Unrecognized format: {format}')

    ro = {r: i for i, r in enumerate(ref_order)}
    intervals = list()
    with open(interval_file, 'rt') as hdl:
        if format == 'SAF':
            next(hdl)
        for line in hdl:
            fields = line.strip().split('\t')
            intervals.append((fields[ix_[0]], int(fields[ix_[1]]), int(fields[ix_[2]]), fields[ix_[3]]))

    intervals.sort(key=lambda iv: interval_key(iv, ro))
    return intervals


class TrackWindow(object):
    """
    Maintains an 'active window' for a given track and allows for internal queries. This is to
    enable multiple overlaps such as

    Read:                <----------------------------------------------
    Track:          [-----------]                              [-------------------]

    """
    def __init__(self, track, refs):
        self.refs = refs
        self.track = track
        self.track_iter = iter(track)
        first_ = next(self.track_iter)
        self.window = [first_]
        self.done_ = False

    def query(self, loc):
        if self.done_ and len(self.window) == 0:
            return []

        self.window = [int for int in self.window if not is_before(int, loc, self.refs)]  # filter
        while len(self.window) == 0 or not is_after(loc, self.window[-1], self.refs):
            try:
                self.window.append(next(self.track_iter))
            except StopIteration:
                self.done_ = True
                break

        return [iv[-1] for iv in self.window if overlaps(iv, loc)]


class TrackOracle(object):
    """
    Manages multiple interval tracks and enables sequential queries to
    genomic positions
    """
    def __init__(self, refnames, tracks=None):
        self.reference_order_ = {k: i for i, k in enumerate(refnames)}
        self.tracks_ = OrderedDict()
        self.offsets_ = OrderedDict()
        if tracks is not None:
            for track in tracks:
                self.add_track(*track)

    def add_track(self, filepath, tag, fmt):
        self.tracks_[tag] = TrackWindow(load_track(filepath, fmt, self.reference_order_), self.reference_order_)

    def query(self, read):
        loc = (read.reference_name, read.reference_start, read.reference_end)
        tag_hits = dict()
        for track in self.tracks_:
            hits = self.tracks_[track].query(loc)
            if hits:
                tag_hits[track] = ','.join(hits)
        return tag_hits


def main(args):
    # load the bam file
    bam = pysam.AlignmentFile(args.bam)
    # marshall tracks
    track_info = [x.split(':') for x in args.track]  # format :: --track /path/to/file.saf:XT:SAF
    tracks = TrackOracle(bam.references, tracks=track_info)
    # open the output bam file with a copied handle
    outbam = pysam.AlignmentFile(args.outbam, template=bam, mode='w')
    for read in bam.fetch():
        if read.is_mapped:
            tags = tracks.query(read)
            for tag, value in tags.items():
                read.set_tag(tag, value)
        outbam.write(read)

    outbam.close()


def get_args():
    parser = ArgumentParser()
    parser.add_argument('bam', help='The input bam to annotate')
    parser.add_argument('outbam', help='The output tagged bam to write')
    parser.add_argument('--track', help='A track to use, format: <filepath>:<tag>:<fmt>; can specify multiple times',
                        action='append')

    return parser.parse_args()


if __name__ == '__main__':
    main(get_args())
