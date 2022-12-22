"""
Simulate RNA reads from paired-tag data

"""
from argparse import ArgumentParser
import numpy as np
import gzip

def get_args():
    parser = ArgumentParser('pairedtag_rna_sim')
    parser.add_argument('read1_fq', help='The read-1 fastq')
    parser.add_argument('read2_fq', help='The read-2 fastq')
    parser.add_argument('gene_count', help='The gene count file')
    parser.add_argument('tx_count', help='The transcript count file')
    parser.add_argument('--cps', help='The number of cells per sample', default=100, type=int)
    parser.add_argument('--upc', help='The number of UMI per cell', default=500, type=int)
    parser.add_argument('--rpu', help='The number of reads per UMI', default=2, type=int)
    parser.add_argument('--fasta', help='The transcript fasta file', default='gencode.v39.transcripts.fa.gz')
    parser.add_argument('--linkers', help='The linker fasta fasta', default='/home/chartl/repos/pipelines/config/linkers.fa')
    parser.add_argument('--sample_bc', help='The sample barcode fasta', default='/home/chartl/repos/pipelines/config/sample_barcode_4bp.fa')
    parser.add_argument('--well_bc', help='The well barcode fasta', default='/home/chartl/repos/pipelines/config/well_barcode_8bp.fa')
    parser.add_argument('--umi_len', help='The umi length', default=10, type=int)
    parser.add_argument('--read_length', help='The lenght of read1', default=50, type=int)
    parser.add_argument('--seed', help='The random seed', type=int, default=1)

    return parser.parse_args()



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


def load_transcripts(tx_fasta):
    gene2tx, tx2gene = dict(), dict()
    for name, seq in read_fasta(tx_fasta):
        if len(seq) < 250:
            continue
        tx, gn, _ = name.split('|', 2)
        tx, gn = tx.split('.')[0], gn.split(".")[0]
        if gn not in gene2tx:
            gene2tx[gn] = list()
        gene2tx[gn].append((tx, seq))
        tx2gene[tx] = gn

    return gene2tx, tx2gene
        
    

M_ = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
      'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
      'N': 'N', 'n': 'n'}

def revC(s):
    return ''.join([M_[x] for x in s])[::-1]


def generate_reads(count_vec, gene2tx, sbarcode, w1bc, w2bc, linkers, rpu, rlen, lumi, rng):
    for cnt, gene in zip(count_vec, gene2tx.keys()):
        if cnt > 0:
           vtx = gene2tx[gene]
           txi = rng.randint(0, len(vtx), cnt)
           txs = [vtx[i] for i in txi]
           for tx in txs:
               umi = ''.join(rng.choice(['A', 'C', 'G', 'T'], lumi))
               q = ''.join(rng.choice(['', 'A', 'C', 'G', 'T'], 3))
               starts = rng.randint(0, len(tx[1]) - rlen, rpu)
               for i, start in enumerate(starts):
                   seq1 = tx[1][start:(start+rlen)]
                   seq2 = umi + w1bc[1] + q + linkers[0][1] + w2bc[1] + linkers[1][1] + 'A' + sbarcode[1] + 'GCGGCCGC'
                   if start % 2 == 1:
                       seq1 = revC(seq1)
                   yield seq1, seq2, i == 0, tx[0], gene


def main(args):
    rng = np.random.RandomState(args.seed)
    gene_transcripts, tx2gene = load_transcripts(args.fasta)
    n_genes = len(gene_transcripts)
    
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


    r1h, r2h, gch, tch = gzip.open(args.read1_fq, 'wt'), gzip.open(args.read2_fq, 'wt'), \
                         gzip.open(args.gene_count, 'wt'), gzip.open(args.tx_count, 'wt') 

    gch.write('id\tsample\twell1\twell2\tgene\tcount\n')
    tch.write('id\tsample\twell1\twell2\tgene\ttranscript\tcount\n')
   
    i_read = 1
    for i_sample, sbc in enumerate(sample_bc):
        i_umi = 1
        i_w1, i_w2 = 0, 0
        cell_counts = rng.multinomial(args.upc, pvals=gene_expr, size=args.cps)
        for i_cell in range(cell_counts.shape[0]):
            read_gen = generate_reads(cell_counts[i_cell,:],
                                      gene_transcripts,
                                      sbc, well_bc[i_w1], well_bc[i_w2], 
                                      linkers, args.rpu, args.read_length, args.umi_len, rng)
            tx_counts, gene_counts = dict(), dict()
            sn, w1n, w2n = sbc[0], well_bc[i_w1][0], well_bc[i_w2][0]
            cid=f'{sn}.{w1n}.{w2n}'
            for r1_seq, r2_seq, is_new_umi, tx, gene in read_gen:
                rname = f'read_{i_read}.u_{i_umi} {sn} {w1n} {w2n} {tx} {gene}'
                r1 = [f'@{rname}', r1_seq, '+', 'I' * len(r1_seq)]
                r2 = [f'@{rname}', r2_seq, '+', 'I' * len(r2_seq)]
                r1h.write('\n'.join(r1) + '\n')
                r2h.write('\n'.join(r2) + '\n')
                if is_new_umi:
                    gene_counts[gene] = 1 + gene_counts.get(gene, 0)
                    tx_counts[tx] = 1 + tx_counts.get(tx, 0)
                    i_umi += 1
                i_read += 1
                if i_read % 1000000 == 0:
                    print(f'Simulated {i_read} reads')

            for gene, count in gene_counts.items():
                gch.write(f'{cid}\t{sn}\t{w1n}\t{w2n}\t{gene}\t{count}\n')

            for tx, count in tx_counts.items():
                txgene = tx2gene[tx]
                txgene = tx2gene[tx]
                txgene = tx2gene[tx]
                txgene = tx2gene[tx]
                tch.write(f'{cid}\t{sn}\t{w1n}\t{w2n}\t{txgene}\t{tx}\t{count}\n')
            i_w2 += 1
            if i_w2 >= len(well_bc):
                i_w1 += 1
                i_w2 = 0

            if i_w1 >= len(well_bc):
                print('Exhausted well barcodes')
                break


if __name__ == "__main__":
    main(get_args())
