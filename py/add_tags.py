"""
Extract barcode information from read names, and place them into
read tags for downstream processing.
"""
from argparse import ArgumentParser
from collections import Counter, OrderedDict
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

        if len(self.window) > 20 or self.done_:
            self.window = [iv for iv in self.window if not is_before(iv, loc, self.refs)]  # filter
        while len(self.window) == 0 or not is_after(self.window[-1], loc, self.refs):
            try:
                self.window.append(next(self.track_iter))
            except StopIteration:
                self.done_ = True
                break

        return [iv[-1] for iv in self.window if overlaps(loc, iv)]


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


def get_args():
    parser = ArgumentParser('add_tags')
    parser.add_argument('bam', help='the input bam file')
    parser.add_argument('out', help='the output bam file')
    parser.add_argument('--library', help='The library name', default='')
    parser.add_argument('--antibody', help='The antibody name', default='')
    parser.add_argument('--track', help='A track to use, format: <filepath>:<tag>:<fmt>; can specify multiple times',
                        action='append')

    return parser.parse_args()


def main(args):
    handle = pysam.AlignmentFile(args.bam)
    sample_ids = Counter(
      (x.query_name.split('|')[2].split(':')[-2],
       x.query_name.split('|')[1].split(':')[-2])
      for x in handle)
    handle.close()

    best = dict()
    for (seq, sam), n in sample_ids.items():
        if sam not in best:
            best[sam] = (seq, n)
        elif n > best[sam][1]:
            best[sam] = (seq, n)

    sample_ids = {k: v[0] for k, v in best.items()}

    handle = pysam.AlignmentFile(args.bam)
    header_dct = handle.header.as_dict()

    rgpfx = args.bam[:-len('.bam')]
    header_dct['RG'] = [
      {'ID': f'{rgpfx}_{i}',
       'PL': 'ILLUMINA',
       'SM': sm,
       'BC': sq, # the BC is not written by pysam even though it's a fine 'alternate' info - add to PU
       'PU': sq} 
      for i, (sm, sq) in enumerate(sample_ids.items())
    ]

    sam2rg = {e['SM']: e['ID'] for e in header_dct['RG']}
    if len(args.track) > 0:
        track_info = [x.split(':') for x in args.track]  # format :: --track /path/to/file.saf:XT:SAF
        tracks = TrackOracle(handle.references, tracks=track_info)
    else:
        tracks = None

    outfile = pysam.AlignmentFile(args.out, mode='wb', header=header_dct)
    for i, read in enumerate(handle):
        outfile.write(transform_read(read, i, sam2rg, args.library, args.antibody, tracks=tracks))
    outfile.close()


def transform_read(read, read_n, sbc_rg_map, library, antibody, tracks=None):
    """
    Given a read of the form
    
    sequencer_read_name|umi:bc1:bc2:sbc:type|umi:bc1:bc2:sbc:type
    
    where the first set of barcodes is parsed, and the second raw,
    extract the UMI and barcode information and place it into
    the expected read tags as follows:
      MI -> umi
      BC -> raw sbc
      CR -> full barcode string
      CB -> parsed BC1BC2
      FC -> the full cell name (lib:anti:sam:bc1:bc2)
      XX -> the read number

    Also add the read to the appropriate read group
 
    Inputs:
    -----------
    read - a pysam AlignedSegment
    read_n - the number of the read
    sbc_rg_map - a dictionary of sample ids to read grop

    Outputs:
    -----------
    transformed read - name reverted to the sequencer name, and tags updated
    """
    rawname, parsed, full = read.query_name.split('|')
    read.set_tag('CR', full)
    read.set_tag('BC', full.split(':')[-2])
    totname_prefix = (library + ':' + antibody).strip(':').lstrip(':')
    psplit = parsed.split(':')
    if totname_prefix != '':
        totname = totname_prefix + f':{psplit[3]}:{psplit[1]}:{psplit[2]}'
    else:
        totname = f'{psplit[3]}:{psplit[1]}:{psplit[2]}'
    read.set_tag('CB', totname)
    read.set_tag('RG', sbc_rg_map[parsed.split(':')[-2]])
    umi = parsed.split(':')[0]
    if umi != '*':
        #if fragmentation happens prior to barcoding, uncomment
        #pos_hash = str(read.reference_start)[-4:]
        #if len(pos_hash) < 4:
        #    pos_hash = '0' * (4 - len(pos_hash)) + pos_hash
        read.set_tag('MI', parsed.split(':')[0]) # + pos_hash)
    read.set_tag('XX', f'{read_n:09d}')
    read.query_name = rawname
    if tracks is not None:
        tags = tracks.query(read)
        for tag, value in tags.items():
            read.set_tag(tag, value)
    return read


if __name__ == '__main__':
    main(get_args())
