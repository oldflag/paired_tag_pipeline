"""
Produce cluster-based QC plots from a PairedTag assay. Be sure to pass in
the `.analysis` h5ad file, which contains merged RNA and DNA counts for
all of the objects.

"""
from argparse import ArgumentParser
import numpy as np
import scanpy as sc
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import cluster_utils as cu


def get_args():
    parser = ArgumentParser()
    parser.add_argument('dna', help='The DNA analysis h5ad file')
    parser.add_argument('rna', help='The RNA analysis h5ad file')
    parser.add_argument('pdf', help='The output QC pdf')
    parser.add_argument('--batch', help='The batch for which to correct', default=None)
    parser.add_argument('--fallback_select', help='Only use count thresholds for selection', action='store_true')
    parser.add_argument('--dna_out', help='base name for the output DNA CSV files', default=None)
    parser.add_argument('--fallback_umi', help='The fallback UMI', default=350, type=int)
    parser.add_argument('--min_rna_umi', type=int, help='Minimum RNA UMI', default=350)
    parser.add_argument('--max_rna_umi', type=int, help='Maximum RNA UMI', default=5000)
    parser.add_argument('--min_dna_umi', type=int, help='Minimum DNA UMI', default=500)
    parser.add_argument('--rna_ngenes', help='Number of genes to use for clustering', default=1500, type=int)
    parser.add_argument('--min_cells_rna', help='Minimum number of cells for RNA', default=500, type=int)

    return parser.parse_args()


def main(args):
    with PdfPages(args.pdf) as pdf:
        cu.PDF = pdf  # set the output
        # rna first
        obj = sc.read_h5ad(args.rna)
        if (obj.obs.antibody_name == 'NA').any():
            dna_fallback = 0
        else:
            dna_fallback = args.fallback_umi
        obj = cu.select_cells(obj, fallback_rna_umi=args.fallback_umi, fallback_dna_umi=dna_fallback, fallback_only=args.fallback_select)
        obj = obj[obj.obs.keep_cell_].copy()
        obj = cu.cluster_pairedtag_rna(obj, min_umi=args.min_rna_umi, n_genes=args.rna_ngenes, max_umi=args.max_rna_umi, harmonize=args.batch,
                                       n_pcs=20, min_cells=args.min_cells_rna)
        # no minimum since select cells did it already
        
        # dna next
        obj = sc.read_h5ad(args.dna)
        obj = obj[np.where(obj.obs.antibody_name != 'NA')[0], :]
        obj = cu.select_cells(obj, fallback_dna_umi=args.fallback_umi, fallback_rna_umi=args.fallback_umi, fallback_only=args.fallback_select)
        obj = obj[obj.obs.keep_cell_].copy()
        obj = cu.cluster_pairedtag_dna(obj, lim_features=1, lim_molecule=args.min_dna_umi, max_molecule=100000, harmonize=args.batch,
                                       n_pcs=10, min_cells_bin=50)

        if args.dna_out is not None:
            for ab in obj.keys():
                of = args.dna_out + '.' + ab + '.csv'
                obj[ab].obs.to_csv(of)


if __name__ == '__main__':
    main(get_args())
