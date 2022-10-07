from argparse import ArgumentParser
from utils import read_fastq
import distance
import numpy as np
import logomaker
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

lseq = 'GTGGCCGATGTTTCGGTGCGAACTCAGACC'
lseq += 'N' * 8
lseq += 'ATCCACGTGCTTGAGAGGCCAGAGCATTCG'

loffset = 18

def get_args():
    p = ArgumentParser()
    p.add_argument('fastq', help='R2 fastq')
    p.add_argument('fasta', help='Linker fasta')
    p.add_argument('logo', help='The logo output file (pdf)')
    return p.parse_args()


def find_best_seq(seq):
    alts = (seq[(loffset+y):(loffset+y+len(lseq))] for y in range(3))
    return min(alts, key=lambda s: distance.hamming(s, lseq))


_lu = {
  'A': np.array([1, 0, 0, 0]),
  'C': np.array([0, 1, 0, 0]),
  'G': np.array([0, 0, 1, 0]),
  'T': np.array([0, 0, 0, 1]),
  'N': np.array([0.25, 0.25, 0.25, 0.25]),
}

def s2c(s):
    return np.vstack([_lu[x] for x in s]) 


def main(args):
    reads = read_fastq(args.fastq)
    PWM = np.zeros((len(lseq), 4), dtype=np.float32)
    with open(args.fasta, 'wt') as out:
        for read in reads:
            out.write('>' + read.name + '\n')
            bs = find_best_seq(read.seq)
            out.write(bs + '\n')
            PWM += s2c(bs)
    PWM = PWM/np.sum(PWM[0,:])
    # transform to information content
    ICM = 2 + np.nansum(PWM * np.log2(PWM), axis=1)
    bitdf = pd.DataFrame(ICM[:, np.newaxis] * PWM)
    bitdf.columns = ('A', 'C', 'G', 'T')
    with PdfPages(args.logo) as pdf:
        logo = logomaker.Logo(df=bitdf,
                              color_scheme='classic',
                              vpad=.1,
                              width=.8)
        logo.style_xticks(rotation=90, fmt='%d', anchor=0)
        pdf.savefig()
    


if __name__ == '__main__':
    main(get_args())
