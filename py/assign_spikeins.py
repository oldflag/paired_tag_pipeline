from argparse import ArgumentParser
from collections import defaultdict
import pysam

def _zcount():
    return {'primary': 0, 'spikein': 0}

def get_args():
    parser = ArgumentParser()
    parser.add_argument('primary_bam')
    parser.add_argument('spike_bam')
    parser.add_argument('filt_primary_bam')
    parser.add_argument('filt_spike_bam')
    parser.add_argument('--info', help='output info file')
    parser.add_argument('--min_qual', type=int, default=30)
    parser.add_argument('--seqtype', help='the sequence type', default='seq')
    parser.add_argument('--dump_counts', help='dump raw counts to this file', default=None)
    parser.add_argument('--duplicate_reads', help='Assign reads to both primary and spikein and do not filter')

    return parser.parse_args()


def iterpairs(readsorted_bam):
    rsb = iter(readsorted_bam)
    for r1, r2 in zip(rsb, rsb):
        yield (r1, r2)

def good_mapping(r1, r2, min_qual):
    either_good = (r1.is_mapped and r1.mapping_quality >= min_qual) or \
                  (r2.is_mapped and r2.mapping_quality >= min_qual)
    both_sanity = True
    if r1.is_mapped and r1.mapping_quality >= min_qual and r2.is_mapped and r2.mapping_quality >= min_qual:
        both_sanity = r1.reference_name == r2.reference_name
    return either_good and both_sanity


def equivalent_mapping(a_r1, a_r2, b_r1, b_r2):
    if good_mapping(a_r1, a_r2, 10) and good_mapping(b_r1, b_r2, 10):
        if (a_r1.is_mapped and a_r2.is_mapped) and (b_r1.is_mapped and b_r2.is_mapped):
            return True
    return False


def mqual(r1, r2):
    q1 = -10 if not r1.is_mapped else r1.mapping_quality
    q2 = -10 if not r2.is_mapped else r2.mapping_quality
    return q1 + q2


def superior_mapping(a_r1, a_r2, b_r1, b_r2):
    return mqual(a_r1, a_r2) > mqual(b_r1, b_r2)
    

def main(args):
    barcode_counts = defaultdict(_zcount)
    pbam = pysam.AlignmentFile(args.primary_bam)
    sbam = pysam.AlignmentFile(args.spike_bam)
    for (r1_prime, r2_prime), (r1_spike, r2_spike) in zip(iterpairs(pbam), iterpairs(sbam)):
        barcode = ':'.join(r1_prime.query_name.split('|')[-1].split(':')[1:3])
        barcode2 = ':'.join(r1_spike.query_name.split('|')[-1].split(':')[1:3])
        assert barcode == barcode2, (barcode, barcode2, r1_prime.query_name, r1_spike.query_name)
        pg, sg = good_mapping(r1_prime, r2_prime, args.min_qual), good_mapping(r1_spike, r2_spike, args.min_qual)
        if pg and not sg:
            barcode_counts[barcode]['primary'] += 1
        elif sg and not pg:
            barcode_counts[barcode]['spikein'] += 1
    pbam.close(), sbam.close()
    barcode_probas = {bc: barcode_counts[bc]['primary']/(0.05+sum(barcode_counts[bc].values())) for bc in barcode_counts}
    pbam = pysam.AlignmentFile(args.primary_bam)
    sbam = pysam.AlignmentFile(args.spike_bam)
    opbam = pysam.AlignmentFile(args.filt_primary_bam, mode='w', template=pbam)
    osbam = pysam.AlignmentFile(args.filt_spike_bam, mode='w', template=sbam)
    for (r1_prime, r2_prime), (r1_spike, r2_spike) in zip(iterpairs(pbam), iterpairs(sbam)):
        barcode = ':'.join(r1_prime.query_name.split('|')[-1].split(':')[1:3])
        if args.duplicate_reads or barcode not in barcode_probas:   # branch 1: keep the reads in both bams
            opbam.write(r1_prime)
            opbam.write(r2_prime)
            osbam.write(r1_spike)
            osbam.write(r2_spike)
        else:   # branch 2: filter the reads
            if barcode_probas.get(barcode, 0) >= 0.7:
                opbam.write(r1_prime)
                opbam.write(r2_prime)
            elif barcode_probas.get(barcode, 1) <= 0.3:
                osbam.write(r1_spike)
                osbam.write(r2_spike)
            else:
                opbam.write(r1_prime)
                opbam.write(r2_prime)
                osbam.write(r1_spike)
                osbam.write(r2_spike)
    opbam.close()
    osbam.close() 

    if args.dump_counts is not None:
        with open(args.dump_counts, 'wt') as out:
            for bc, tab in barcode_counts.items():
                out.write(bc + '\t' + str(tab['primary']) + '\t' + str(tab['spikein']) + \
                          '\t%.1f%%\t%.1f%%\t%.3f\n' % (
                           (100 * tab['primary'])/(tab['primary'] + tab['spikein']),
                           (100 * tab['spikein']/(tab['primary'] + tab['spikein'])),
                           barcode_probas.get(bc, -1)))

    barcode_tots = {bc: sum(cn.values()) for bc, cn in barcode_counts.items()}

    ohdl = open(args.info, 'wt')
    spname = args.primary_bam.split('__')
    lib, sam, ab = spname[:3]
    ohdl.write('library,sample,antibody,library_type,count_threshold,tot_bcs,primary_barcodes,primary_reads,spikein_barcodes,spikein_reads,ambiguous_barcodes\n')
    for count_thr in (0, 50, 100, 200, 500):
        spike_reads, prime_reads = 0, 0
        tot_bcs, spike_bcs, prime_bcs, ambig_bcs = 0, 0, 0, 0
        for bc, tot_cnt in barcode_tots.items():
            if tot_cnt < count_thr: 
                continue
            tot_bcs += 1
            spike_reads += barcode_counts[bc]['spikein']
            prime_reads += barcode_counts[bc]['primary']
            if barcode_probas[bc] > 0.8:
                prime_bcs += 1
            elif barcode_probas[bc] < 0.2:
                spike_bcs += 1
            else:
                ambig_bcs += 1
        ohdl.write(f'{lib},{sam},{ab},{args.seqtype},{count_thr},{tot_bcs},{prime_bcs},{prime_reads},{spike_bcs},{spike_reads},{ambig_bcs}\n')
    ohdl.close()


if __name__ == '__main__':
    main(get_args())
