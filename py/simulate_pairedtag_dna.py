"""
Simulate RNA reads from paired-tag data

"""
from argparse import ArgumentParser
import numpy as np
import gzip
import os

def get_args():
    parser = ArgumentParser('pairedtag_dna_sim')
    parser.add_argument('read1_fq', help='The read-1 fastq')
    parser.add_argument('read2_fq', help='The read-2 fastq')
    parser.add_argument('peak_count', help='The peak count file')
    parser.add_argument('--npeaks', help='The number of peaks to simulate', default=50000)
    parser.add_argument('--cps', help='The number of cells per sample', default=100, type=int)
    parser.add_argument('--upc', help='The number of UMI per cell', default=500, type=int)
    parser.add_argument('--rpu', help='The number of reads per UMI', default=2, type=int)
    parser.add_argument('--fasta', help='The genome fasta file; must be indexed', default='/home/share/storages/2T/genome/human/GRCh38.primary_assembly.genome.fa')
    parser.add_argument('--linkers', help='The linker fasta fasta', default='/home/chartl/repos/pipelines/config/linkers.fa')
    parser.add_argument('--sample_bc', help='The sample barcode fasta', default='/home/chartl/repos/pipelines/config/sample_barcode_4bp.fa')
    parser.add_argument('--well_bc', help='The well barcode fasta', default='/home/chartl/repos/pipelines/config/well_barcode_8bp.fa')
    parser.add_argument('--umi_len', help='The umi length', default=10, type=int)
    parser.add_argument('--read_length', help='The lenght of read1', default=50, type=int)
    parser.add_argument('--seed', help='The random seed', type=int, default=1)

    args = parser.parse_args()
    if not os.path.exists(args.fasta + '.fai'):
        raise ValueError('FASTA file must have a .fai file')


    return args 



def read_fasta(tx_fasta):
    if tx_fasta[-2:] == 'gz':
        hdl = gzip.open(tx_fasta, 'rt')
    else:
        hdl = open(tx_fasta, 'rt')

    name, seq = None, None
    for line in hdl:
        if line[0] == '>':
            if name:
                yield name, seq
            name, seq = line[1:-1], ''
        else:
            seq += line[:-1]
    yield name, seq


def get_peaks(genome_fasta, n, gen):
    fai = genome_fasta + '.fai'
    contigs = list()
    for entry in (x.strip().split() for x in open(fai, 'rt')):
        if entry[0][:3] == 'chr':
            contigs.append((entry[0], int(entry[1])))

    w = np.array([x[1] for x in contigs])
    w = w/np.sum(w)

    ppc = gen.multinomial(int(1.075*n), w)
    peaks = list()
    for (contig, size), npeaks in zip(contigs, ppc):
        start, end = 1000000, size - 1000000
        peak_starts = gen.randint(start, end, npeaks)
        peak_ends = peak_starts + gen.randint(95,175, npeaks) * (1 + gen.binomial(3, 0.15, npeaks))
        peaks.extend([(contig, ps, pe) for ps, pe in zip(peak_starts, peak_ends)]) 

    with open('sim_peaks.bed', 'wt') as out:
        for i, (c, s, e) in enumerate(peaks):
            out.write(f'{c}\t{s}\t{e}\tpeak_{1+i}\t.\t.\n')

    os.system(f'bedtools getfasta -name+ -bed sim_peaks.bed -fi {genome_fasta} -fo sim_peaks.fa')
    fasta_entries = read_fasta('sim_peaks.fa')
    fasta_entries = [x for x in fasta_entries if 'N' not in x[1]]  # remove 'N' sequences
    if len(fasta_entries) > n:
        return fasta_entries[:n]

    print('Warning:: after filtering N bases, only %d peaks remain' % len(fasta_entries))
    return fasta_entries
        
    

M_ = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
      'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
      'N': 'N', 'n': 'n'}

def revC(s):
    return ''.join([M_[x] for x in s])[::-1]


def generate_reads(count_vec, peak_list, sbarcode, w1bc, w2bc, linkers, rpu, rlen, lumi, rng):
    for cnt, peak in zip(count_vec, peak_list):
        if cnt > 0:
           for _ in range(cnt):
               umi = ''.join(rng.choice(['A', 'C', 'G', 'T'], lumi))
               q = ''.join(rng.choice(['', 'A', 'C', 'G', 'T'], 3))
               start = rng.randint(0, len(peak[1]) - rlen, 1)[0]
               seq1 = peak[1][start:(start+rlen)]
               seq2 = umi + w1bc[1] + q + linkers[0][1] + w2bc[1] + linkers[1][1] + 'A' + sbarcode[1] + 'GCGGCCGC'
               if start % 2 == 1:
                   seq1 = revC(seq1)
               for rn in range(rpu):
                   yield seq1, seq2, rn == 0, peak[0]
    

def main(args):
    rng = np.random.RandomState(args.seed)
    peaks = get_peaks(args.fasta, args.npeaks, rng)
    n_genes = len(peaks)
    
    sample_bc = list(read_fasta(args.sample_bc))
    well_bc = list(read_fasta(args.well_bc))
    linkers = list(read_fasta(args.linkers))
    nsam = len(sample_bc)
    gene_expr = rng.dirichlet(alpha=np.ones(n_genes))

    tu = args.upc * args.cps * len(sample_bc)
    tr = args.rpu * tu
    sstr = (
        'Will simulate:\n' +
       f'   {len(sample_bc)} samples\n' +
       f'   {args.cps} cells/sample\n' +
       f'   {len(sample_bc)*args.cps} total cells\n' +
       f'   {args.upc} UMI per cell\n' +
       f'   {args.rpu} reads per umi\n' +
       f'   {tu} / {tr}  total UMI / reads\n\n')

    print(sstr)


    r1h, r2h, pch = gzip.open(args.read1_fq, 'wt'), gzip.open(args.read2_fq, 'wt'), \
                    gzip.open(args.peak_count, 'wt')

    pch.write('id\tsample\twell1\twell2\tpeak\tcount\n')
   
    i_read = 1
    for i_sample, sbc in enumerate(sample_bc):
        i_umi = 1
        i_w1, i_w2 = 0, 0
        cell_counts = rng.multinomial(args.upc, pvals=gene_expr, size=args.cps)
        for i_cell in range(cell_counts.shape[0]):
            read_gen = generate_reads(cell_counts[i_cell,:],
                                      peaks,
                                      sbc, well_bc[i_w1], well_bc[i_w2], 
                                      linkers, args.rpu, args.read_length, args.umi_len, rng)
            peak_counts = dict()
            sn, w1n, w2n = sbc[0], well_bc[i_w1][0], well_bc[i_w2][0]
            cid=f'{sn}.{w1n}.{w2n}'
            for r1_seq, r2_seq, is_new_umi, gene in read_gen:
                rname = f'read_{i_read}.u_{i_umi} {sn} {w1n} {w2n} {gene}'
                r1 = [f'@{rname}', r1_seq, '+', 'I' * len(r1_seq)]
                r2 = [f'@{rname}', r2_seq, '+', 'I' * len(r2_seq)]
                r1h.write('\n'.join(r1) + '\n')
                r2h.write('\n'.join(r2) + '\n')
                if is_new_umi:
                    peak_counts[gene] = 1 + peak_counts.get(gene, 0)
                    i_umi += 1
                i_read += 1
                if i_read % 1000000 == 0:
                    print(f'Simulated {i_read} reads')

            for gene, count in peak_counts.items():
                pch.write(f'{cid}\t{sn}\t{w1n}\t{w2n}\t{gene}\t{count}\n')

            i_w2 += 1
            if i_w2 >= len(well_bc):
                i_w1 += 1
                i_w2 = 0

            if i_w1 >= len(well_bc):
                print('Exhausted well barcodes')
                break


if __name__ == "__main__":
    main(get_args())
