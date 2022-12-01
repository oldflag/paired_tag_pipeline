from argparse import ArgumentParser
import pysam
import numpy as np
from multiprocessing import Pool


def get_args():
    parser = ArgumentParser()
    parser.add_argument('bam', help='the bam file')
    parser.add_argument('sai', help='the annotation (gene) sai file')
    parser.add_argument('out', help='the output score file')
    parser.add_argument('--ncore', help='the number of cores', type=int, default=1)

    return parser.parse_args()


def parse_sai(saif):
    with open(saif, 'rt') as inf:
        next(inf)
        for line in inf:
            fields = line.strip().split('\t')
            yield fields[0], fields[1], int(fields[2]), int(fields[3]), fields[4]


def calc_tsse_(bampath, gene_locus):
    name, contig, gen_start, gen_end, strand = gene_locus
    if strand == '+':
        rgn_s = max(gen_start - 1000, 1)
        rgn_e = gen_start + 1000
    else:
        rgn_s = max(1, gen_end - 1000)
        rgn_e = gen_end + 1000

    bam = pysam.AlignmentFile(bampath)
    cvg = dict()
    bcp = set()
    for read in bam.fetch(contig, rgn_s, rgn_e):
        ucp = str(read.reference_start) + read.get_tag('CB') + (read.get_tag('MI')[0] if read.has_tag('MN') else 'NN')
        if ucp not in bcp:
            antibody = read.get_tag('RG').split('__')[-2]
            if antibody not in cvg:
                cvg[antibody] = np.zeros(shape=(2000,), dtype=int)
            bcp.add(ucp)
            ref_start, ref_end = read.reference_start, read.reference_end
            array_start = max(0, ref_start - rgn_s)
            array_end = min(2000, array_start + (ref_end - ref_start))
            cvg[antibody][array_start:array_end] = 1 + cvg[antibody][array_start:array_end]
    # generate the bins
    ret = dict()
    for antibody in cvg:
        bins = np.array([np.sum(cvg[antibody][x:(x+100)]) for x in range(1900)])
        #th = max(bins[:4].mean(), bins[-4:].mean())
        th = (bins[0] + bins[-1])/2
        if th == 0:
            if all(bins == 0):
                th = 1
            else:
                th = np.min(bins[bins > 0])
        ret[antibody] = (name, np.mean(bins), bins/th)
    return ret


def calc_tsse(args):
    return calc_tsse_(*args)


def main(args):
    gene_locs = parse_sai(args.sai)
    if args.ncore > 1:
        pool = Pool(args.ncore)
        gene_info = pool.imap(calc_tsse, ((args.bam, loc) for loc in gene_locs))
    else:
        gene_info = map(calc_tsse, ((args.bam, loc) for loc in gene_locs))

    with open(args.out, 'wt') as out:
        fieldnames = ['gene', 'antibody', 'depth', 'tsse'] + ['bin_%d' % (i+1) for i in range(1900)]
        out.write('\t'.join(fieldnames) + '\n')
        for info in gene_info:
            for antibody, (gene, mean_dp, bins) in info.items():
                escore = np.max(bins)
                fields = [gene, antibody, str(mean_dp), str(escore)] + [str(x) for x in bins]
                out.write('\t'.join(fields) + '\n')
    
    pool.close()


if __name__ == '__main__':
    main(get_args())
