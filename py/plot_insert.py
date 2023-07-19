from argparse import ArgumentParser
import pandas as pd
import numpy as np
from plotnine import ggplot, geom_point, geom_smooth, aes, geom_bar
from plotnine import facet_grid, save_as_pdf_pages, facet_wrap



def get_args():
    parser = ArgumentParser()
    parser.add_argument('fragfile')
    parser.add_argument('out_pdf')
    parser.add_argument('--primary', help='Name of the primary species', default='human')

    return parser.parse_args()


def main(ffile, opdf, primary_name):
    dat_insert = pd.read_csv(ffile, sep='\t', header=0, names=['bam', 'insert', 'count'])
    if '__' in dat_insert.bam[0]:
        dat_insert['library'] = dat_insert.bam.astype(str).apply(lambda s: s.split('__')[0])
        dat_insert['antibody'] = dat_insert.bam.astype(str).apply(lambda s: s.split('__')[2])
        dat_insert['sampleid'] = dat_insert.bam.astype(str).apply(lambda s: s.split('__')[1][:8])
        other_name = 'mouse' if primary_name == 'human' else 'human'
        dat_insert['species'] = dat_insert.bam.astype(str).apply(lambda s: other_name if 'spike' in s else primary_name)
        dat_insert = dat_insert[dat_insert.antibody != 'UNK']
    else:
        # bulk
        dat_insert['antibody'] = dat_insert.bam.astype(str).apply(lambda s: s.split('_')[0])
        dat_insert['sampleid'] = 'bulk'
        dat_insert['library'] = 'bulk'
        dat_insert['species'] = primary_name

    insert_libs = dat_insert.groupby(['antibody', 'library', 'species', 'insert'])['count'].apply(sum).reset_index()
    insert_tots = insert_libs.groupby(['antibody', 'library', 'species'])['count'].apply(sum).reset_index()
    insert_tots.columns = ['antibody', 'library', 'species', 'total_count']

    pbar = ggplot(insert_tots) + geom_bar(aes(x='antibody', y='total_count', fill='species'),
                                  stat='identity', position='dodge')
    insert_mg = insert_libs.merge(insert_tots, on=['antibody', 'species', 'library'])
    insert_mg['normalized_count'] = insert_mg['count']/insert_mg['total_count']
    insert_mg['abspecies'] = insert_mg['antibody'] + ':' + insert_mg['species']
    bscatter = ggplot(insert_mg[insert_mg['insert'] < 1000]) + \
               geom_point(aes(x='insert', y='normalized_count', color='antibody'),
                          size=0.1) + \
               facet_grid("antibody~species", scales='free')

    iprim = insert_mg[insert_mg['species'] == primary_name]
    pscatter = ggplot(iprim[iprim['insert'] < 1000]) + \
               geom_point(aes(x='insert', y='normalized_count', color='antibody'),
                          size=0.1) + \
               facet_wrap('~ antibody')

    
    isec = insert_mg[insert_mg['species'] != primary_name]
    iscatter = ggplot(isec[(isec['insert'] < 1000) & (isec['normalized_count'] < np.percentile(isec['normalized_count'], 99))]) + \
               geom_point(aes(x='insert', y='normalized_count', color='antibody'),
                          size=0.1) + \
               facet_wrap('~ antibody')

    save_as_pdf_pages([pbar, bscatter, pscatter, iscatter], filename=opdf)


if __name__ == '__main__':
    args = get_args()
    main(args.fragfile, args.out_pdf, args.primary)
