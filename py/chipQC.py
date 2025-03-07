"""
UMI-aware single-cell ChIPQC

Given:
  - An aligned, sorted .bam file
  - A peak .bed file

Produces the following metrics:
  - total reads
  - avg reads / cell
  - median reads / cell
  - std reads / cell
  - aligned reads
  - avg aligned / cell
  - median aligned /cell
  - std aligned / cell
  - # of called peaks
  - peak size
  - peak size sd
  - mean peak cvg
  - reads in peak
  - reads in peak / cell

Using the EE and EP tags:
  - reads in enhancer
  - reads in enhancer / cell
  - reads in promoter
  - reads in promoter / cell

Using the GN (gene name) tags:
  - reads in genebody
  - reads in genebody / cell

"""

from argparse import ArgumentParser
from collections import Counter, defaultdict
from functools import reduce
from multiprocessing import Pool
import itertools
import pysam
import sys


STANDARD_CONTIGS = {'%d' % x for x in range(1, 26)} | {'X', 'Y'}
STANDARD_CONTIGS |= {'chr%s' % u for u in STANDARD_CONTIGS}

def safediv(a, b):
    if b == 0:
        return 0
    else:
        return a/b


def get_args():
    parser = ArgumentParser('chipQC')
    parser.add_argument('bam', help='The bam file')
    parser.add_argument('peaks', help='The peak file')
    parser.add_argument('out', help='The output stats file')
    parser.add_argument('--verbose', help='Be loud', action='store_true')
    parser.add_argument('--sample_out', help='The sample statistics output file', default=None)

    return parser.parse_args()


def getUMI(read):
    return read.get_tag('MI') if read.has_tag('MI') else None


def group_umi(reads):
    rbg = defaultdict(list)
    for read in reads:
        if read.has_tag('RG'):
            grp = (read.get_tag('RG'), read.get_tag('CB'), getUMI(read))
        else:
            lib = read.get_tag('CB').split(':')[0]
            grp = (lib + '__UNK__UNK__UNK', read.get_tag('CB'), read.get_tag('CR').split(':')[0]) 
        rbg[grp].append(read)
    return rbg


def iter_umi(samfile, verb=False):
    """
    Iterate over UMIs within a bam file. This is done by:

    1) Aggregating all reads at a position
    2) Splitting reads from (1) by the barcode and UMI
    3) Yielding each group from (2)

    """
    read = next(samfile)
    while read.is_unmapped:
        read = next(samfile)
    if read.is_unmapped:
        raise ValueError('No reads mapped or SAM file out of order')
    loc = (read.reference_name, read.reference_start)
    stash = [read]
    for i, read in enumerate(samfile):
        if read.is_unmapped:
            if len(stash) > 0:
                yield loc, group_umi(stash)
                break
        nloc = (read.reference_name, read.reference_start)
        if loc[1] != nloc[1] or loc[0] != nloc[0]:
            yield loc, group_umi(stash)
            stash, loc = [read], nloc
        else:
            stash.append(read)
        if verb and ((i+1) % 1000000) == 0:
            sys.stderr.write('Processed %d reads\n' % (i+1))
    else:  # no unmapped reads
        if len(stash) > 0:
            yield loc, group_umi(stash)


def cstats():
    return {
      'filtered_MQ30': {'in_peak': {'reads': 0, 'umi': 0},
                        'off_target': {'reads': 0, 'umi': 0},
                        'in_enhancer': {'reads': 0, 'umi': 0},
                        'in_promoter': {'reads': 0, 'umi': 0},
                        'in_genebody': {'reads': 0, 'umi': 0}
                       },
      'retained': {'in_peak': {'reads': 0, 'umi': 0},
                   'off_target': {'reads': 0, 'umi': 0},
                   'in_enhancer': {'reads': 0, 'umi': 0},
                   'in_promoter': {'reads': 0, 'umi': 0},
                   'in_genebody': {'reads': 0, 'umi': 0}
                  },
      'unmapped': {'reads': 0},
      'no_UMI': {'reads': 0}
    }


def read_peaks(peakfile, refs):
    with open(peakfile, 'rt') as hdl:
        next(hdl)   # header
        peaks = [(x[1], int(x[2]), int(x[3])) for x in (l.strip().split('\t') for l in hdl)]  # name,chr,start,stop,strand
        peaks.sort(key=lambda x: (refs.index(x[0]), x[1], x[2]))
        for peak in peaks:
            yield peak


def ends_before(read, loc, refs):
    if read.reference_name != loc[0]:
        return refs.index(read.reference_name) < refs.index(loc[0])
    else:
        return read.reference_end < loc[1]


def starts_after(read, loc, refs):
    if read.reference_name != loc[0]:
        return refs.index(read.reference_name) > refs.index(loc[0])
    else:
        return read.reference_start > loc[2]


def loc_starts_after(rloc, loc, refs):
    if rloc[0] != loc[0]:
        return refs.index(rloc[0]) > refs.index(loc[0])
    else:
        return rloc[1] > loc[2]

def overlaps(read, loc, refs):
    return not (ends_before(read, loc, refs) or starts_after(read, loc, refs))


def next_peak(peak, rloc, peak_iter, refs):
    while loc_starts_after(rloc, peak, refs):
        try:
            peak = next(peak_iter) 
        except StopIteration:
            peak = (refs[-1], 10 ** 10, 10 ** 10 + 1)
    return peak


def div_(a, b):
    if b < 1e-9:
        return 0
    return a/b


def compute_sample_stats(cell_stats, a2s_map):
    sample_stats = dict()
    for fc in cell_stats:
        if '*' not in fc:
            library, antibody, assay, w1, w2 = fc.split(':') 
        else:
            library, antibody, assay = fc.split(':')[:3]
            continue   # can't assign the cell, don't count it
        # re-key
        sample_id = a2s_map.get(assay, assay)
        sample = library + ':' + sample_id + ':' + antibody
        good_reads = cell_stats[fc]['retained']['in_peak']['reads'] + \
                     cell_stats[fc]['retained']['off_target']['reads']
        good_umi = cell_stats[fc]['retained']['in_peak']['umi'] + \
                   cell_stats[fc]['retained']['off_target']['umi']
        good_on_target_reads = cell_stats[fc]['retained']['in_peak']['reads']
        good_on_target_umi = cell_stats[fc]['retained']['in_peak']['umi'] 
        filtered_reads = cell_stats[fc]['filtered_MQ30']['in_peak']['reads'] + \
                         cell_stats[fc]['filtered_MQ30']['off_target']['reads']
        filtered_umi = cell_stats[fc]['filtered_MQ30']['in_peak']['umi'] + \
                         cell_stats[fc]['filtered_MQ30']['off_target']['umi']
        unmapped_reads = cell_stats[fc]['unmapped']['reads']
        if sample not in sample_stats:
            sample_stats[sample] = {
                'total_reads': list(),
                'HQ_reads': list(),
                'total_umi': list(),
                'HQ_umi': list(),
                'umi_in_peak': list(),
                'reads_in_peak': list(),
                'estimated_unmapped_umi': list(),
                'reads_in_enhancer': list(),
                'reads_in_promoter': list(),
                'reads_in_genebody': list(),
                'umi_in_enhancer': list(),
                'umi_in_promoter': list(),
                'umi_in_genebody': list(),
                'reads_in_CRE': list(),
                'umi_in_CRE': list()
            }
        sample_stats[sample]['total_reads'].append(good_reads + filtered_reads + unmapped_reads)
        sample_stats[sample]['total_umi'].append(good_umi + filtered_umi)
        sample_stats[sample]['estimated_unmapped_umi'].append((good_umi/good_reads) * unmapped_reads if good_reads > 0 else -1)
        sample_stats[sample]['HQ_umi'].append(good_umi)
        sample_stats[sample]['HQ_reads'].append(good_reads)
        sample_stats[sample]['umi_in_peak'].append(good_on_target_umi)
        sample_stats[sample]['reads_in_peak'].append(good_on_target_reads)
        sample_stats[sample]['reads_in_enhancer'].append(cell_stats[fc]['retained']['in_enhancer']['reads'])
        sample_stats[sample]['umi_in_enhancer'].append(cell_stats[fc]['retained']['in_enhancer']['umi'])
        sample_stats[sample]['reads_in_promoter'].append(cell_stats[fc]['retained']['in_promoter']['reads'])
        sample_stats[sample]['umi_in_promoter'].append(cell_stats[fc]['retained']['in_promoter']['umi'])
        sample_stats[sample]['reads_in_genebody'].append(cell_stats[fc]['retained']['in_genebody']['reads'])
        sample_stats[sample]['umi_in_genebody'].append(cell_stats[fc]['retained']['in_genebody']['umi'])
        sample_stats[sample]['reads_in_CRE'].append(sample_stats[sample]['reads_in_enhancer'][-1] + \
                                                    sample_stats[sample]['reads_in_promoter'][-1])
        sample_stats[sample]['umi_in_CRE'].append(sample_stats[sample]['umi_in_enhancer'][-1] + \
                                                  sample_stats[sample]['umi_in_promoter'][-1])

    # define some further statistics
    DERIVED_STATISTICS = {
        'read_peak_rate': (('reads_in_peak', 'HQ_reads'), lambda x: safediv(x[0],x[1])),
        'umi_peak_rate': (('umi_in_peak', 'HQ_umi'), lambda x: safediv(x[0], x[1])),
        'read_enhancer_rate': (('reads_in_enhancer', 'HQ_reads'), lambda x: safediv(x[0], x[1])),
        'umi_enhancer_rate': (('umi_in_enhancer', 'HQ_umi'), lambda x: safediv(x[0], x[1])),
        'read_promoter_rate': (('reads_in_promoter', 'HQ_reads'), lambda x: safediv(x[0], x[1])),
        'umi_promoter_rate': (('umi_in_promoter', 'HQ_umi'), lambda x: safediv(x[0], x[1])),
        'read_CRE_rate': (('reads_in_CRE', 'HQ_reads'), lambda x: safediv(x[0], x[1])),
        'umi_CRE_rate': (('umi_in_CRE', 'HQ_reads'), lambda x: safediv(x[0], x[1])),
        'read_genic_rate': (('reads_in_genebody', 'HQ_reads'), lambda x: safediv(x[0], x[1])),
        'umi_genic_rate': (('umi_in_genebody', 'HQ_umi'), lambda x: safediv(x[0], x[1])),
        'enhancer_promoter_ratio': (('umi_in_enhancer', 'umi_in_promoter'), lambda x: safediv(x[0], x[1]))
    } 
    sample_vals = dict()
    for sample in sample_stats:
        sample_vals[sample] = dict()
        umi_vec = sample_stats[sample]['total_umi']
        for key in sample_stats[sample]:
            cell_values = sample_stats[sample][key]
            total = sum(sample_stats[sample][key]) 
            nobs = len(cell_values)
            cell_values = [x for x, c in zip(cell_values, umi_vec) if c >= 750]
            filtered_total = sum(cell_values)
            m = div_(sum(cell_values), len(cell_values))
            vsorted = sorted(cell_values)
            i_25, i_50, i_75 = int(0.25*len(cell_values)), int(0.5*len(cell_values)), int(0.75*len(cell_values))
            if len(cell_values) > 0:
                sample_vals[sample][key] = {
                   'total': total,
                   'n_cell': nobs,
                   'n_well_covered': len(cell_values),
                   'well_covered_total': filtered_total,
                   'min': min(cell_values),
                   'Q25': vsorted[i_25],
                   'median': vsorted[i_50],
                   'Q75': vsorted[i_75],
                   'max': max(cell_values),
                   'mean': m,
                   'std': sum(((x-m)**2 for x in cell_values))/len(cell_values)
                }
            else:
                sample_vals[sample][key] = {
                   'n_cell': nobs,
                   'total': total,
                   'n_well_covered': 0,
                   'well_covered_total': 0,
                   'min': 'NA',
                   'Q25': 'NA',
                   'median': 'NA',
                   'Q75': 'NA',
                   'max': 'NA',
                   'mean': 'NA',
                   'std': 'NA'
                }

        for statname, (keylst, statfx) in DERIVED_STATISTICS.items():
            key_lists = [sample_stats[sample][k] for k in keylst]
            cell_values = [statfx(x) for x in zip(*key_lists)]
            total = statfx([sample_vals[sample][k]['total'] for k in keylst])
            wc_total = statfx([sample_vals[sample][k]['well_covered_total'] for k in keylst])
            nobs = len(cell_values)
            cell_values = [x for x, c in zip(cell_values, umi_vec) if c >= 750]
            filtered_total = sum(cell_values)
            m = div_(filtered_total, len(cell_values))
            vsorted = sorted(cell_values)
            i_25, i_50, i_75 = int(0.25*len(cell_values)), int(0.5*len(cell_values)), int(0.75*len(cell_values))
            if len(cell_values) > 0:
                sample_vals[sample][statname] = {
                    'total': total,
                    'n_cell': nobs,
                    'n_well_covered': len(cell_values),
                    'well_covered_total': filtered_total,
                    'min': min(cell_values),
                    'Q25': vsorted[i_25],
                    'median': vsorted[i_50],
                    'Q75': vsorted[i_75],
                    'max': max(cell_values),
                    'mean': m,
                    'std': sum(((x-m)**2 for x in cell_values))/len(cell_values)
                }
            else:
                sample_vals[sample][statname] = {
                   'n_cell': nobs,
                   'total': total,
                   'n_well_covered': 0,
                   'well_covered_total': 0,
                   'min': 'NA',
                   'Q25': 'NA',
                   'median': 'NA',
                   'Q75': 'NA',
                   'max': 'NA',
                   'mean': 'NA',
                   'std': 'NA'
                }

    return sample_vals
   


def main(args):
    """
    Compute ChIP-seq statistics on a per-cell and aggregate level. For each read:
      (i) The UI is extracted from the MI tag and deduplicated at that location
      (ii) The cell is extracted from the CB tag
         + Update total reads, total umi
         + (if aligned, standard contig) update aligned reads, aligned umi
         + (if in peak) update peak reads, peak umi
    """
    bam_hdl = pysam.AlignmentFile(args.bam)
    # the metadata has been encoded in the read group lines; pull out the assay -> sample map
    assay2sample = dict()
    for rg_ifo in bam_hdl.header['RG']:
        rg = rg_ifo['ID']
        seq, antibody, assay, sample, _ = rg.split('__')
        assay2sample[assay] = sample

    peak_iter = read_peaks(args.peaks, bam_hdl.references)
    peak = next(peak_iter)
    stat_counts = defaultdict(cstats)
    for loc, umi_group in iter_umi(bam_hdl, args.verbose):
        ppeak = peak
        peak = next_peak(peak, loc, peak_iter, bam_hdl.references)
        for (libsample, cell, umi), reads in umi_group.items():
            fcell = cell
            if umi is None:
                stat_counts[fcell]['no_UMI']['reads'] += len(reads)
                continue
            u = False
            k2 = 'in_peak' if overlaps(reads[0], peak, bam_hdl.references) else 'off_target'
            for read in reads:
                k1 = 'retained' if read.mapq >= 30 else 'filtered_MQ30'
                stat_counts[fcell][k1][k2]['reads'] += 1
                if read.has_tag('EP'):
                    stat_counts[fcell][k1]['in_promoter']['reads'] += 1
                    stat_counts[fcell][k1]['in_promoter']['umi'] += (1 if u is False else 0)
                elif read.has_tag('EE'):
                    stat_counts[fcell][k1]['in_enhancer']['reads'] += 1
                    stat_counts[fcell][k1]['in_enhancer']['umi'] += (1 if u is False else 0)
                if read.has_tag('GN'):
                    stat_counts[fcell][k1]['in_genebody']['reads'] += 1
                    stat_counts[fcell][k1]['in_genebody']['umi'] += (1 if u is False else 0)
                if k1 == 'retained' and u is False:
                    stat_counts[fcell][k1][k2]['umi'] += 1
                    u = True
            if u is False:
                stat_counts[fcell]['filtered_MQ30'][k2]['umi'] += 1
    for read in bam_hdl:   # now just unmapped
        cell = read.get_tag('CB')
        stat_counts[fcell]['unmapped']['reads'] += 1
                 

    with open(args.out, 'wt') as out:
        out.write('sample_id\tlibrary\tsample\tantibody\tassay\tcell\tfilter\ttarget\ttype\tcount\n')
        for fullcell in stat_counts:
            if '*' not in fullcell:
                library, antibody, assay_id, w1, w2 = fullcell.split(':')
                cell = w1 + ":" + w2
            else:
                library, antibody, assay_id = fullcell.split(':')[:3]
                cell = '*'
            assay2sample.get(assay_id, assay_id)
            for k1 in ('retained', 'filtered_MQ30'):
                for k2 in ('in_peak', 'off_target', 'in_enhancer', 'in_promoter', 'in_genebody'):
                    for k3 in ('reads', 'umi'):
                        cct = stat_counts[fullcell][k1][k2][k3]
                        out.write(f'{assay_id}\t{library}\t{sample}\t{antibody}\t{assay_id}\t{cell}\t{k1}\t{k2}\t{k3}\t{cct}\n')
            n_unmapped = stat_counts[fullcell]['unmapped']['reads']
            out.write(f'{assay_id}\t{library}\t{sample}\t{antibody}\t{assay_id}\t{cell}\tunmapped\tNA\treads\t{n_unmapped}\n')

   
    if args.sample_out is not None:
        # compute per-sample statistics
        with open(args.sample_out, 'wt') as out:
            out.write('sample_id\tlibrary\tsample\tantibody\tfeature\tmetric\tvalue\n')
            for sample_id, sample_stats in compute_sample_stats(stat_counts, assay2sample).items():
                library, sample, antibody = sample_id.split(':')[:3]
                for k1, dct in sample_stats.items():
                    for k2, v in dct.items():
                        out.write(f'{sample_id}\t{library}\t{sample}\t{antibody}\t{k1}\t{k2}\t{v}\n')


if __name__ == '__main__':
    main(get_args())
