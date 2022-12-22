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

def safeget(dct, itm):
    return dct.get(itm, 'NA')


def get_args():
    parser = ArgumentParser('count2h5')
    parser.add_argument('counts', help='Output 3-column counts file (e.g. from UMI count)')
    parser.add_argument('sample_digest', help='The sample digest CSV file containing sequence ID and library information, will be added to `.obs`')
    parser.add_argument('h5', help='The output h5ad file')
    parser.add_argument('--split_feature', help='In the case of multiple features, redistribute single counts proportionally', action='store_true')

    return parser.parse_args()


def get_sid(bc_string):
    return bc_string.split('.')[-1]


def main(args):
    row_map, row_ind = OrderedDict(), list()
    col_map, col_ind = OrderedDict(), list()
    cell_counts, cell_features = dict(), dict()
    data = list()
    with xopen(args.counts, 'rt') as hdl:
        next(hdl)  # header
        for feature, cell, count_str in (x.strip().split('\t') for x in hdl):
            if '*' in cell:
                continue
            if cell not in row_map:
                row_map[cell] = len(row_map)
                cell_counts[cell], cell_features[cell] = 0, 0
            cell_counts[cell] += int(count_str)
            cell_features[cell] += 1
            if ',' not in feature:
                if feature not in col_map:
                    col_map[feature] = len(col_map)
                data.append(int(count_str))
                row_ind.append(row_map[cell])
                col_ind.append(col_map[feature])
            else:
                subfeatures = feature.split(',')
                n_feats = len(subfeatures)
                if args.split_feature:
                    newcount = float(count_str)/n_feats
                else:
                    newcount = int(count_str)
                for sf in subfeatures:
                    if sf not in col_map:
                        col_map[sf] = len(col_map)
                    data.append(newcount)
                    row_ind.append(row_map[cell])
                    col_ind.append(col_map[sf])

    data_type = float if args.split_feature else int
               
    data, row_ind, col_ind = np.array(data, dtype=data_type), np.array(row_ind, dtype=int), np.array(col_ind, dtype=int)
    sample_info = pd.read_csv(args.sample_digest)
    sample_info.loc[:, 'sample_key'] = sample_info.assay_info + '_' + sample_info.sequence_id

    X = scipy.sparse.csr_matrix((data, (row_ind, col_ind)))
    col_meta = pd.DataFrame({'feature_name': list(col_map.keys())})
    row_meta = pd.DataFrame({'atom_id': list(row_map.keys()),
                             'feature_count': [cell_features[x] for x in row_map],
                             'molecule_count': [cell_counts[x] for x in row_map]})
    print(row_meta.head().to_string())
    row_meta.loc[:, 'sequence_id'] = row_meta.atom_id.map(lambda x: x.split(':')[0])
    row_meta.loc[:, 'antibody_name'] = row_meta.atom_id.map(lambda x: x.split(':')[1])
    row_meta.loc[:, 'sample_id'] = row_meta.atom_id.map(lambda x: x.split(':')[2])
    row_meta.loc[:, 'well1_id'] = row_meta.atom_id.map(lambda x: x.split(':')[3])
    row_meta.loc[:, 'well2_id'] = row_meta.atom_id.map(lambda x: x.split(':')[4])
    row_meta.loc[:, 'sample_key'] = row_meta.sample_id + '_' + row_meta.sequence_id
    print(row_meta.head().to_string())
    # add in metadata
    for cname in sample_info.columns:
        if cname in {'sequence_id', 'antibody_name', 'sample_id', 'sample_key'}:
            continue
        mapping = dict(zip(sample_info.sample_key, sample_info.loc[:, cname]))
        row_meta.loc[:, cname] = np.array(row_meta.sample_key.map(lambda x: safeget(mapping, x)).values, dtype=str)
    row_meta.loc[:, 'cell_id'] = row_meta.sample_id + ':' + row_meta.antibody_name + ':' + row_meta.lysis_id + ':' + row_meta.well1_id + ':' + row_meta.well2_id
    print(row_meta.head().to_string())
    for cn in row_meta.columns:
        print('%s -> %s' % (cn, row_meta.loc[:, cn].dtype.__repr__()))
    anndat = AnnData(X, row_meta, col_meta)
    anndat.write_h5ad(args.h5, compression='gzip')


if __name__ == '__main__':
    main(get_args())
