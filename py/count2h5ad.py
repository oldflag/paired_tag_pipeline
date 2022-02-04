"""
Transforms a 3-column count file of the form
  <gene>   <cell>   <#>
  <gene>   <cell>   <#>

into a compressed sparse row matrix, and inserts the matrix,
rownames, and columnnames into an AnnData object, and saves
it as an h5ad file.

"""
from argparse import ArgumentParser
from anndata import AnnData
import scipy
import pandas as pd
import numpy as np
from utils import xopen
from collections import OrderedDict

def get_args():
    parser = ArgumentParser('count2h5')
    parser.add_argument('counts', help='Output 3-column counts file (e.g. from UMI count)')
    parser.add_argument('h5', help='The output h5ad file')

    return parser.parse_args()


def get_sid(bc_string):
    i = bc_string.index('SM')
    return bc_string[i:]


def main(args):
    row_map, row_ind = OrderedDict(), list()
    col_map, col_ind = OrderedDict(), list()
    data = list()
    with xopen(args.counts, 'rt') as hdl:
        next(hdl)  # header
        for feature, cell, count_str in (x.strip().split('\t') for x in hdl):
            if cell not in row_map:
                row_map[cell] = len(row_map)
            if feature not in col_map:
                col_map[feature] = len(col_map)
            data.append(int(count_str))
            row_ind.append(row_map[cell])
            col_ind.append(col_map[feature])
    data, row_ind, col_ind = np.array(data, dtype=int), np.array(row_ind, dtype=int), np.array(col_ind, dtype=int)

    X = scipy.sparse.csr_matrix((data, (row_ind, col_ind)))
    col_meta = pd.DataFrame({'feature_name': list(col_map.keys())})
    row_meta = pd.DataFrame({'cell_id': list(row_map.keys()),
                             'sample_id': [get_sid(x) for x in row_map]})
    print(row_meta.head().to_string())

    anndat = AnnData(X, row_meta, col_meta)
    anndat.write_h5ad(args.h5, compression='gzip')


if __name__ == '__main__':
    main(get_args())
