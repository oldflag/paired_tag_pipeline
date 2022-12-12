"""
Given a csv file with cluster IDs and tag information, split into cluster-specific bam files
"""
from argparse import ArgumentParser
import pysam
import pandas as pd
import numpy as np

def get_args():
    parser = ArgumentParser()
    parser.add_argument('bam', help='The bam file')
    parser.add_argument('csv', help='The csv file (metadata)')
    parser.add_argument('--cluster_col', help='use this column for cluster ids')
    parser.add_argument('--out_dir', help='the output directory', default='.')

    return parser.parse_args()


def main(args):
    # load the atom -> csv map
    dat = pd.read_csv(args.csv)

    base_name = args.bam[:-4]  # drop ".bam"

    out_handles = dict()
    clusters = dat.loc[:, args.cluster_col].astype(str).unique()
    na_ix = np.where(pd.isna(dat.antibody))[0]
    in_handle = pysam.AlignmentFile(args.bam)
    hdr = pysam.AlignmentFile(args.bam, 'rb')
    dat.loc[:, 'antibody'][na_ix] = 'NA'
    antibodies = dat.antibody.astype(str).unique()
    for antibody in antibodies:
        dat_sub = dat[dat.antibody == antibody]
        atom_to_clust = dict(zip(dat_sub.atom_id, dat_sub.loc[:, args.cluster_col].astype(str)))
        print('Processing %s (%d atoms)' % (antibody, len(atom_to_clust)))
        for cname in clusters:
            bf = f'{args.out_dir}/{base_name}_{args.cluster_col}_{cname}.{antibody}.bam'
            hdl = pysam.AlignmentFile(bf, header=in_handle.header, mode='wb')
            dat_clst = dat_sub[dat_sub.loc[:, args.cluster_col] == cname]
            for atm in dat_clst.atom_id:
                out_handles[atm] = hdl
    
    for i, read in enumerate(in_handle.fetch()):
        atom = read.get_tag('CB', None)
        if atom in out_handles:
            out_handles[atom].write(read)
        if i % 20000000 == 0:
            print('.... ' + str(i))

    print('Processed %d reads on %s' % (i, antibody))
    in_handle.close()
            
    for k, v in out_handles.items():
        v.close()
        
    in_handle.close()

if __name__ == '__main__':
    main(get_args())
