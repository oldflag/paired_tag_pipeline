"""
Produce cluster-based QC plots from a PairedTag assay. Be sure to pass in
the `.analysis` h5ad file, which contains merged RNA and DNA counts for
all of the objects.

"""
from argparse import ArgumentParser
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

    return parser.parse_args()


def main(args):
    with PdfPages(args.pdf) as pdf:
        cu.PDF = pdf  # set the output
        # rna first
        obj = sc.read_h5ad(args.rna)
        obj = cu.select_cells(obj, fallback_rna_umi=300, fallback_dna_umi=300)
        obj = obj[obj.obs.keep_cell_].copy()
        obj = cu.cluster_pairedtag_rna(obj, min_umi=1, n_genes=2000, max_umi=3000,
                                       harmonize=args.batch, n_pcs=15)  
        # no minimum since select cells did it already
        
        # dna next
        obj = sc.read_h5ad(args.dna)
        obj = cu.select_cells(obj, fallback_dna_umi=300, fallback_rna_umi=300)
        obj = obj[obj.obs.keep_cell_].copy()
        obj = cu.cluster_pairedtag_dna(obj, lim_features=1, lim_molecule=1, max_molecule=4000,
                                       harmonize=args.batch)




if __name__ == '__main__':
    main(get_args())
