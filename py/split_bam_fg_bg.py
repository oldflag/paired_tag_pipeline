import pysam
from argparse import ArgumentParser
import numpy as np
from collections import defaultdict
from multiprocessing import Pool

HELP_ = """\
Takes in any number of single-cell Cut&Tag bam files, and writes reads to
separate "foreground" and "background" bam files based on the barcode count.
Input bams should be sorted by position; output bams are not required to be sorted.

For each input .bam file:
  + The "CB" (cell barcode) and "MI" (umi) tags of the reads are used to count total
    UMI per barcode
  + The TLIM is set to the 99%ile of the UMI for the top 3000 cell barcodes
  + The LLIM is set to 10% of the TLIM
  + "Foreground" barcodes are defined as those with total counts >= LLIM
  + "Background" barcodes are defined as those with total counts < LLIM
  + On a second pass through the .bam file, all reads corresponding to a "Foreground" CB
    are written to the foreground .bam; otherwise to the background .bam
"""

def get_args():
    parser = ArgumentParser(HELP_)
    parser.add_argument('--fg_out', help='The output foreground bam file')
    parser.add_argument('--bg_out', help='The output background bam file')
    parser.add_argument('--threads', help='Number of threads', default=8, type=int)
    parser.add_argument('bam', help='Input bam file(s)', nargs='+')

    return parser.parse_args()


def get_cb(read):
    return read.get_tag('CB') if read.has_tag('CB') else ':'.join(read.query_name.split('|')[1].split(':')[1:-1])


def get_mi(read):
    return read.get_tag('MI') if read.has_tag('MI') else read.query_name.split('|')[1].split(':')[0]


def count_umis_in_bam(bam_file):
    local_umi_counts = defaultdict(lambda: defaultdict(int))
    
    with pysam.AlignmentFile(bam_file, 'rb') as infile:
        for read in infile:
            if read.is_unmapped:  # Skip unmapped reads
                continue
            rb = get_cb(read)
            cb = (bam_file + '_' + rb)
            mi = get_mi(read)
            contig = read.reference_name
            pos = read.reference_start
            if cb and mi:
                unique_umi = (contig, pos, mi)
                local_umi_counts[cb][unique_umi] += 1
                
    return {k: {k2: v2 for k2, v2 in v.items()} for k, v in local_umi_counts.items()}


def reduce_counts(counts_list):
    merged_counts = dict()
    for dict_ in counts_list:
        for k, v in dict_.items():
            merged_counts[k] = v
    
    return merged_counts


def main(args):
    # Initialize
    barcode_umi_counts = defaultdict(lambda: defaultdict(int))
    tot_sams = len(args.bam)
    exp_cells = int(tot_sams * (4000 / 12))  # conservative estimate of cell recovery

    # Count UMIs for each cell barcode
    with Pool(processes=args.threads) as pool:
        counts_list = pool.map(count_umis_in_bam, args.bam)

    barcode_umi_counts = reduce_counts(counts_list)
    # Get top 3000 cell barcodes based on total UMI count
    barcode_umi_sums = {cb: len(umis) for cb, umis in barcode_umi_counts.items()}
    top_barcodes = sorted(barcode_umi_sums, key=barcode_umi_sums.get, reverse=True)[:exp_cells]
    top_umi_counts = [barcode_umi_sums[cb] for cb in top_barcodes]
    tumi = np.array([t for t in barcode_umi_sums.values()])

    # Calculate TLIM and LLIM
    TLIM = np.percentile(np.array(top_umi_counts), 99)
    LLIM = 0.1 * TLIM

    print('Cutting at %d UMI (%d barcodes)' % (LLIM, np.sum(tumi >= LLIM)))

    # Initialize BAM output files
    fg_out = pysam.AlignmentFile(args.fg_out, 'wb', template=pysam.AlignmentFile(args.bam[0]))
    bg_out = pysam.AlignmentFile(args.bg_out, 'wb', template=pysam.AlignmentFile(args.bam[0]))

    # Second pass to write to foreground and background BAM files
    for bam_file in args.bam:
        with pysam.AlignmentFile(bam_file, 'rb') as infile:
            for read in infile:
                cb = bam_file + '_' + get_cb(read)
                if cb and barcode_umi_sums.get(cb, 0) >= LLIM:
                    fg_out.write(read)
                elif cb and barcode_umi_sums.get(cb,0) <= LLIM * 0.2:
                    bg_out.write(read)

    # Close BAM files
    fg_out.close()
    bg_out.close()


if __name__ == '__main__':
    main(get_args())
