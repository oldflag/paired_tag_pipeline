"""
Extract barcode information from read names, and place them into
read tags for downstream processing.
"""
from argparse import ArgumentParser
from collections import Counter
import pysam


def get_args():
    parser = ArgumentParser('add_tags')
    parser.add_argument('bam', help='the input bam file')
    parser.add_argument('out', help='the output bam file')
    parser.add_argument('--library', help='The library name', default='')
    parser.add_argument('--antibody', help='The antibody name', default='')

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

    outfile = pysam.AlignmentFile(args.out, mode='wb', header=header_dct)
    for i, read in enumerate(handle):
        outfile.write(transform_read(read, i, sam2rg, args.library, args.antibody))
    outfile.close()


def transform_read(read, read_n, sbc_rg_map, library, antibody):
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
    return read


if __name__ == '__main__':
    main(get_args())
