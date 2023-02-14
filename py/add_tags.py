"""
Extract barcode information from read names, and place them into
read tags for downstream processing.
"""
from argparse import ArgumentParser
from collections import Counter, OrderedDict
from csv import DictReader
import pysam


def interval_key(loc1, ctg_odr):
    return 10**9 * ctg_odr[loc1[0]] + loc1[1] + 10**(-6) * (loc1[2] - loc1[1])


def overlaps(loc1, loc2, ctg_odr=None):
    return (loc1[0] == loc2[0]) and (loc2[2] >= loc1[1]) and (loc2[1] <= loc1[2])


def is_before(loc1, loc2, ctg_odr):
    #print((loc1, loc2))
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
            if fields[ix_[0]] in ro:
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
        if loc is None or loc[1] == -1 or loc[0] not in self.refs:
            return []
        if self.done_ and len(self.window) == 0:
            return []
        if len(self.window) > 20 or self.done_:
            self.window = [iv for iv in self.window if iv[1] != -1 and not is_before(iv, loc, self.refs)]  # filter
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
                tag_hits[track] = ','.join(sorted(set(hits)))
        return tag_hits


def get_args():
    parser = ArgumentParser('add_tags')
    parser.add_argument('bam', help='the input bam file')
    parser.add_argument('out', help='the output bam file')
    parser.add_argument('--library', help='The library name', default='')
    parser.add_argument('--antibody', help='The antibody name', default='')
    parser.add_argument('--assay_id', help='The assay id', default=None)
    parser.add_argument('--sample_id', help='The sample id', default=None)
    parser.add_argument('--track', help='A track to use, format: <filepath>:<tag>:<fmt>; can specify multiple times',
                        action='append')
    parser.add_argument('--reachtools', help='Flag that the bam was produced by reachtools', action='store_true')
    parser.add_argument('--sample_digest', help='The sample digest file - only needed with --reachtools', default=None)
    parser.add_argument('--library_digest', help='The library digest file - only needed with --reachtools', default=None)

    return parser.parse_args()


def parse_samples_pp(args):
    """\
    Parse samples from a pipelines-produced bam file
    """
    if args.assay_id is None:
        handle = pysam.AlignmentFile(args.bam)
        sample_ids = Counter(
          (x.query_name.split('|')[2].split(':')[-2],
           x.query_name.split('|')[1].split(':')[-2])
          for x in handle if x.query_name[-1] != '*')
        handle.close()

        best = dict()
        for (seq, sam), n in sample_ids.items():
            if sam not in best:
                best[sam] = (seq, n)
            elif n > best[sam][1]:
                best[sam] = (seq, n)

        sample_ids = {k: v[0] for k, v in best.items()}
    else:
        sample_ids = {'USER_SPEC': args.assay_id}
    return sample_ids, None


def parse_samples_rt(args):
    """\
    Parse samples from a reachtools-produced bam file and aux files
    """
    assert args.sample_digest is not None
    assert args.library_digest is not None
    # read the library digest
    lib_records = [x for x in DictReader(open(args.library_digest), delimiter=',')]
    # get the sequence id from the bam file
    bam2sid = {rec['bam'].split('/')[-1]: rec['sequence_id'] for rec in lib_records}
    seq_id = bam2sid[args.bam.split('/')[-1]]
    # for reachtools we want a map from the well to the assay id
    assay2well, well2info = dict(), dict()
    for rec in DictReader(open(args.sample_digest), delimiter=','):
        if rec['sequence_id'] == seq_id:
            wid = rec['dna_well']
            wpad = ('0' if len(wid) < 2 else '') + wid
            assay2well[rec['assay_id']] = wpad
            well2info[wpad] = rec
    return assay2well, well2info


def parse_samples(args):
    if args.reachtools:
        return parse_samples_rt(args)
    return parse_samples_pp(args)
    

def main(args):
    print(args)
    sample_ids, rt_assay_info = parse_samples(args)  # rt_assay_info only for ReachTools
    handle = pysam.AlignmentFile(args.bam)
    print('Loading header')
    header_dct = handle.header.as_dict()

    suffix = None
    if args.assay_id is not None:
        suffix = args.assay_id
    if args.sample_id is not None:
        if suffix:
            suffix += '__' + args.sample_id
            sample_ids = {args.sample_id: args.assay_id}
        else:
            suffix = '__' + args.sample_id

    rgpfx = args.bam[:-len('.bam')] if suffix is None else args.library + '__' + args.antibody + '__' + suffix
    header_dct['RG'] = [
      {'ID': f'{rgpfx}__{i}',
       'PL': 'ILLUMINA',
       'SM': sm,
       'BC': sq, # the BC is not written by pysam even though it's a fine 'alternate' info - add to PU
       'PU': sq} 
      for i, (sm, sq) in enumerate(sample_ids.items())
    ]
    print('Header written')

    sam2rg = {e['SM']: e['ID'] for e in header_dct['RG']} if args.assay_id is None else {args.assay_id: f'{rgpfx}_0'}
    if args.track is not None and len(args.track) > 0:
        track_info = [x.split(':') for x in args.track]  # format :: --track /path/to/file.saf:XT:SAF
        tracks = TrackOracle(handle.references, tracks=track_info)
    else:
        tracks = None

    outfile = pysam.AlignmentFile(args.out, mode='wb', header=header_dct)
    transformer = transform_read_rt if args.reachtools else transform_read
    transformer_args = {'assay_info': rt_assay_info} if args.reachtools else dict()
    for i, read in enumerate(handle):
        outfile.write(transformer(read, i, sam2rg, args.library, args.antibody, tracks=tracks, **transformer_args))
    outfile.close()


def transform_read_rt(read, read_n, sbc_rg_map, library, antibody, assay_info, tracks=None):
    """
    Given a read of the form

    (seq:...:prefix):WIDX1:WIDX2:SIDX:UMI

    extract the UMI and barcode information and place it into the expected
    read tags as follows:
      MI -> UMI
      BC -> assay_id
      CR -> full barcode string
      CB -> the full cell name (lib:anti:sam:bc1:bc2)
      AB -> the antibody
      XX -> the read number

    And add the read to the appropriate read group. In addition, annotate
    any overlapping tracks.

    
    Inputs:
    -----------
    read - a pysam AlignedSegment
    read_n - the number of the read
    sbc_rg_map - a dictionary of (assay) ids to read grop
    library - the library to use [ignored]
    antibody - the name of the antibody [ignored]
    assay_info - a dictionary of sample well -> dict() containing the correct information
    tracks - A TrackWindow object holding genomic tracks

    Outputs:
    -----------
    transformed read - name reverted to the sequencer name, and tags updated
    """
    split_fields = read.query_name.split(':')
    bc2, bc1, sm_well, umi = split_fields[-4:]
    antibody, library = assay_info[sm_well]['antibody_name'], assay_info[sm_well]['library_id']
    assay_id, sample_id = assay_info[sm_well]['assay_id'], assay_info[sm_well]['sample_id']
    read.set_tag('BC', '%s:%s:%s:%s' % (umi, bc1, bc2, sm_well))
    read.set_tag('CR', '%s:%s:%s:%s' % (umi, bc1, bc2, assay_id))
    read.set_tag('CB', '%s:%s:%s:%s:%s' % (library, antibody, assay_id, bc1, bc2))
    read.set_tag('MI', umi)
    read.set_tag('AB', antibody)
    read.set_tag('RG', sbc_rg_map[assay_id])
    read.set_tag('XX', f'{read_n:09d}')
    read.query_name = ':'.join(split_fields[:-4])
    if tracks is not None:
        tags = tracks.query(read)
        for tag, value in tags.items():
            read.set_tag(tag, value)
    return read


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
      CB -> the full cell name (lib:anti:sam:bc1:bc2)
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
    try:
        read.set_tag('BC', full.split(':')[-2])
    except:
        pass
    totname_prefix = (library + ':' + antibody).strip(':').lstrip(':')
    psplit = parsed.split(':')
    if totname_prefix != '' and len(psplit) > 2:
        totname = totname_prefix + f':{psplit[3]}:{psplit[1]}:{psplit[2]}'
    elif len(psplit) > 2:
        totname = f'{psplit[3]}:{psplit[1]}:{psplit[2]}'
    elif totname_prefix != '':
        totname = totname_prefix + ':*'
    else:
        totname = '*'
    read.set_tag('CB', totname)
    try:
        read.set_tag('RG', sbc_rg_map[parsed.split(':')[-2]])
        umi = parsed.split(':')[0]
    except:
        umi = '*'
    if umi != '*':
        #if fragmentation happens prior to barcoding, uncomment
        #pos_hash = str(read.reference_start)[-4:]
        #if len(pos_hash) < 4:
        #    pos_hash = '0' * (4 - len(pos_hash)) + pos_hash
        read.set_tag('MI', umi)
    read.set_tag('XX', f'{read_n:09d}')
    #read.query_name = rawname
    if tracks is not None:
        tags = tracks.query(read)
        for tag, value in tags.items():
            read.set_tag(tag, value)
    return read


if __name__ == '__main__':
    main(get_args())
