import seaborn as sbn
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import numpy as np
import scanpy as sc
from argparse import ArgumentParser

def get_args():
    parser = ArgumentParser()
    parser.add_argument('rna_gene_read')
    parser.add_argument('rna_gene_umi')
    parser.add_argument('dna_bin_read')
    parser.add_argument('dna_bin_umi')
    parser.add_argument('outpdf')
    return parser.parse_args()

def main(args):
    # only using metadata - so read the h5ad and extract just `.obs`
    dat_rna_read = sc.read_h5ad(args.rna_gene_read).obs
    dat_rna_umi = sc.read_h5ad(args.rna_gene_umi).obs
    dat_dna_read = sc.read_h5ad(args.dna_bin_read).obs
    dat_dna_umi = sc.read_h5ad(args.dna_bin_umi).obs
    dat_rna_read.rename({'feature_count': 'rna_ngenes_read',
                         'molecule_count': 'rna_reads'}, axis=1, inplace=True)
    dat_rna_umi.rename({'feature_count': 'rna_ngenes_umi',
                        'molecule_count': 'rna_umis'}, axis=1, inplace=True)
    dat_dna_read.rename({'feature_count': 'dna_nbins_read',
                         'molecule_count': 'dna_reads'}, axis=1, inplace=True)
    dat_dna_umi.rename({'feature_count': 'dna_nbins_umi',
                        'molecule_count': 'dna_umis'}, axis=1, inplace=True)

    dat_mg = dat_rna_read.merge(dat_rna_umi[['cell_id', 'rna_ngenes_umi', 'rna_umis']], on='cell_id')
    dat_mg = dat_mg.merge(dat_dna_read[['cell_id', 'dna_nbins_read', 'dna_reads']], on='cell_id')
    dat_mg = dat_mg.merge(dat_dna_umi[['cell_id', 'dna_nbins_umi', 'dna_umis']], on='cell_id')

    # update the RNA and DNA objects to reflect the merged data
    new_cols = ['dna_nbins_read', 'dna_reads', 'dna_nbins_umi',
                'dna_umis', 'rna_ngenes_read', 'rna_reads', 'rna_ngenes_umi',
                'rna_umis']
    dat_rna_umi = sc.read_h5ad(args.rna_gene_umi)
    dat_rna_umi = dat_rna_umi[dat_rna_umi.obs.cell_id.isin(dat_mg.cell_id)].copy()
    rna_mg_df = dat_rna_umi.obs
    rna_mg_df = rna_mg_df.merge(dat_mg[['cell_id'] + new_cols], on='cell_id', how='left')
    for cn in new_cols:
        dat_rna_umi.obs.loc[:, cn] = rna_mg_df.loc[:, cn].values
    rna_umi_fn = args.rna_gene_umi[:-5] + '.analysis.h5ad'
    dat_rna_umi.write_h5ad(rna_umi_fn)

    del dat_rna_umi, rna_mg_df

    dat_dna_umi = sc.read_h5ad(args.dna_bin_umi)
    dat_dna_umi = dat_dna_umi[dat_dna_umi.obs.cell_id.isin(dat_mg.cell_id)].copy()
    dna_mg_df = dat_dna_umi.obs
    dna_mg_df = dna_mg_df.merge(dat_mg[['cell_id'] + new_cols], on='cell_id', how='left')
    for cn in new_cols:
        print(dna_mg_df.loc[:, cn].head())
        dat_dna_umi.obs.loc[:, cn] = dna_mg_df.loc[:, cn].values
        print(dat_dna_umi.obs.loc[:, cn].head())
    dna_umi_fn = args.dna_bin_umi[:-5] + '.analysis.h5ad'
    dat_dna_umi.write_h5ad(dna_umi_fn)

    del dna_mg_df, dat_dna_umi

    library_sizes = dat_mg.groupby(['lysis_id'])[['dna_reads', 'dna_umis', 'rna_reads', 'rna_umis']].sum().reset_index()

    dnar_map = dict(zip(library_sizes.lysis_id, library_sizes.dna_reads))
    dnau_map = dict(zip(library_sizes.lysis_id, library_sizes.dna_umis))
    rnar_map = dict(zip(library_sizes.lysis_id, library_sizes.rna_reads))
    rnau_map = dict(zip(library_sizes.lysis_id, library_sizes.rna_umis))

    dat_mg.loc[:, 'rna_reads_norm'] = dat_mg.rna_reads / dat_mg.lysis_id.map(rnar_map.__getitem__).astype(np.float32)
    dat_mg.loc[:, 'rna_umis_norm'] = dat_mg.rna_umis / dat_mg.lysis_id.map(rnau_map.__getitem__).astype(np.float32)
    dat_mg.loc[:, 'dna_reads_norm'] = dat_mg.rna_umis / dat_mg.lysis_id.map(dnar_map.__getitem__).astype(np.float32)
    dat_mg.loc[:, 'dna_umis_norm'] = dat_mg.dna_umis / dat_mg.lysis_id.map(dnau_map.__getitem__).astype(np.float32)

    dat = dat_mg

    def yxline(min_=None, max_=None):
        ax = plt.gca()
        lims = [
            min_ or np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
            max_ or np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
        ]

        # now plot both limits against eachother
        ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
        ax.set_aspect('equal')
        ax.set_xlim(lims)
        ax.set_ylim(lims)
    with PdfPages(args.outpdf) as pdf: 
        plt.figure(figsize=(12,8))
        for lib in dat.lysis_id.unique():
            ddf = dat[dat.lysis_id == lib]
            plt.scatter(ddf.rna_reads, ddf.rna_umis, marker='.', label=lib, rasterized=True)
        yxline();
        plt.legend(bbox_to_anchor=(1,1))
        plt.tight_layout()
        plt.xlabel('RNA Reads (each cell)')
        plt.ylabel('RNA UMI (each cell)');
        
        plt.gca().set_rasterized(True)
        #plt.savefig(pdf,format='pdf',dpi=250)
        plt.close()
        plt.figure(figsize=(12,8))
        for lib in dat.lysis_id.unique():
            ddf = dat[dat.lysis_id == lib]
            plt.scatter(ddf.rna_reads, ddf.rna_umis, marker='.', label=lib)
        yxline(min_=1)
        plt.xscale('log'); plt.yscale('log');
        plt.legend(bbox_to_anchor=(1,1))
        plt.tight_layout()
        plt.xlabel('RNA Reads (each cell)')
        plt.ylabel('RNA UMI (each cell)'); 
        plt.gca().set_rasterized(True)
        #plt.savefig(pdf,format='pdf',dpi=250)
        plt.close()


        plt.figure(figsize=(12,8))
        for lib in dat.lysis_id.unique():
            ddf = dat[dat.lysis_id == lib]
            plt.scatter(ddf.dna_reads, ddf.dna_umis, marker='.', label=lib)
        yxline();
        plt.legend(bbox_to_anchor=(1,1))
        plt.tight_layout()
        plt.xlabel('DNA Reads (each cell)')
        plt.ylabel('DNA UMI (each cell)');
        
        plt.gca().set_rasterized(True)
        #plt.savefig(pdf,format='pdf',dpi=250)
        plt.close()

        plt.figure(figsize=(12,8))
        for lib in dat.lysis_id.unique():
            ddf = dat[dat.lysis_id == lib]
            plt.scatter(ddf.dna_reads, ddf.dna_umis, marker='.', label=lib)
        yxline(min_=1)
        plt.xscale('log'); plt.yscale('log');
        plt.legend(bbox_to_anchor=(1,1))
        plt.tight_layout()
        plt.xlabel('RNA Reads (each cell)')
        plt.ylabel('RNA UMI (each cell)');
        plt.gca().set_rasterized(True);
        #plt.savefig(pdf,format='pdf',dpi=250)
        plt.close()

        dat.loc[:, 'dna_complexity'] = dat.dna_umis / dat.dna_reads
        dat.loc[:, 'rna_complexity'] = dat.rna_umis / dat.rna_reads
        plt.figure(figsize=(8,6))
        
        for lib in dat.lysis_id.unique():
            ddf = dat[dat.lysis_id == lib]
            ddf = ddf[(ddf.dna_umis >= 500) & (ddf.rna_umis >= 500)]
            if ddf.shape[0] > 0:
                plt.scatter(ddf.rna_complexity, ddf.dna_complexity, marker='.', label=lib, rasterized=True)
        
        plt.legend(bbox_to_anchor=(1,1))
        plt.tight_layout()
        plt.xlabel('RNA Complexity (UMI per read; per cell)')
        plt.ylabel('DNA Complexity (UMI per read; per cell)');
        plt.title('Library Complexity per cell\n(500 UMI minimum per cell)');yxline();
        plt.tight_layout();plt.gca().set_rasterized(True)
        #plt.savefig(pdf,format='pdf',dpi=250)
        plt.close()

        for lib in dat.lysis_id.unique():
            dsub = dat[dat.lysis_id == lib]
            dsub.to_csv('%s.csv' % lib)

            fig, axs = plt.subplots(3,4,constrained_layout=True,figsize=(12,8))
            for i, sample_id in enumerate(dsub.sample_id.unique()):
                ddf = dsub[dsub.sample_id == sample_id]
                ddf = ddf[ddf.rna_reads > 0]
                ddf = ddf[ddf.dna_reads > 0]
                if ddf.shape[0] > 0:
                    col = i % 4
                    row = min(2, int(i/4))
                    ax = axs[row, col]
                    ax.scatter(ddf.rna_reads, ddf.dna_reads, marker='.', label=sample_id, s=0.5)
                    ax.hlines(y=300, xmin=1, xmax=ddf.rna_reads.max(), linestyles='dashed', colors='red', linewidth=0.5)
                    ax.vlines(x=300, ymin=1, ymax=ddf.dna_reads.max(), linestyles='dashed', colors='red', linewidth=0.5)
                    ax.set_xscale('log');ax.set_yscale('log');
                    ax.set_xlabel('RNA Reads (each cell)');ax.set_ylabel('DNA reads (each cell)')
                    ax.set_title(sample_id, size=8)
                    ax.set_rasterized(True)
            fig.suptitle(lib, fontsize=12)
            plt.gca().set_rasterized(True);plt.savefig(pdf, format='pdf', dpi=250);plt.close()
            

            fig, axs = plt.subplots(3,4,constrained_layout=True,figsize=(12,8))
            for i, sample_id in enumerate(dsub.sample_id.unique()):
                ddf = dsub[dsub.sample_id == sample_id]
                ddf = ddf[ddf.rna_umis > 0]
                ddf = ddf[ddf.dna_umis > 0]
                if ddf.shape[0] > 0:
                    col = i % 4
                    row = min(2, int(i/4))
                    ax = axs[row, col]
                    ax.scatter(ddf.rna_umis, ddf.dna_umis, marker='.', label=sample_id, s=0.5)
                    ax.hlines(y=300, xmin=1, xmax=ddf.rna_umis.max(), linestyles='dashed', colors='red', linewidth=0.5)
                    ax.vlines(x=300, ymin=1, ymax=ddf.dna_umis.max(), linestyles='dashed', colors='red', linewidth=0.5)
                    ax.set_xscale('log');ax.set_yscale('log');
                    ax.set_xlabel('RNA UMI (each cell)');ax.set_ylabel('DNA UMI (each cell)')
                    ax.set_title(sample_id, size=8)
                    ax.set_rasterized(True)
            fig.suptitle(lib, fontsize=12)
            plt.gca().set_rasterized(True);plt.savefig(pdf, format='pdf', dpi=250);plt.close()
        

        plt.figure(figsize=(12,12))
        for lib in dat.lysis_id.unique():
            ddf = dat[dat.lysis_id == lib]
            plt.scatter(ddf.rna_umis_norm, ddf.dna_umis_norm, marker='.', label=lib)
        yxline();
        plt.legend(bbox_to_anchor=(1,1))
        plt.tight_layout()
        plt.xlabel('RNA UMIs (normalized)')
        plt.ylabel('DNA UMIs (normalized)');
        yxline(min_=5e-8);
        plt.xscale('log');plt.yscale('log')
        
        
        plt.gca().set_rasterized(True);plt.savefig(pdf,format='pdf',dpi=250);plt.close(); plt.figure(figsize=(12,12))
        lowlab = 'low'
        for lib in dat.lysis_id.unique():
            ddf = dat[dat.lysis_id == lib]
            ddf1 = ddf[(ddf.dna_umis >= 300) & (ddf.rna_umis >= 300)]
            ddf2 = ddf[(ddf.dna_umis < 300)  | (ddf.rna_umis < 300)]
            plt.scatter(ddf2.rna_umis_norm, ddf2.dna_umis_norm, marker='.', label=lowlab, color='grey')
            plt.scatter(ddf1.rna_umis_norm, ddf1.dna_umis_norm, marker='.', label=lib)
            lowlab = '_nolegend_'
        yxline();
        plt.legend(bbox_to_anchor=(1,1))
        plt.tight_layout()
        plt.xlabel('RNA UMIs (normalized)')
        plt.ylabel('DNA UMIs (normalized)');
        yxline(min_=5e-8);
        plt.xscale('log');plt.yscale('log')
        plt.gca().set_rasterized(True)
        #plt.savefig(pdf,format='pdf',dpi=250)
        plt.close()
        
        covered_both = dat[(dat.rna_umis >= 500) & (dat.dna_umis >=500)].groupby('lysis_id')['rna_umis'].count().reset_index()
        covered_both.columns = ['library', 'n_cells']
        covered_both.loc[:, 'assay'] = 'dna+rna'
        covered_rna = dat[(dat.rna_umis >=500)].groupby('lysis_id')['rna_umis'].count().reset_index()
        covered_rna.columns = ['library', 'n_cells']
        covered_rna.loc[:, 'assay'] = 'rna'
        covered_dna = dat[(dat.dna_umis >= 500)].groupby('lysis_id')['dna_umis'].count().reset_index()
        covered_dna.columns = ['library', 'n_cells']
        covered_dna.loc[:, 'assay'] = 'dna'
        
        ddf = pd.concat([covered_rna, covered_dna, covered_both])
        plt.figure(figsize=(12,8))
        sbn.barplot(x='library', y='n_cells', hue='assay', data=ddf)
        plt.xticks(rotation=90);plt.ylabel('Number of covered cells');plt.xlabel('Library')
        plt.tight_layout();plt.gca().set_rasterized(True)
        #plt.savefig(pdf, format='pdf', dpi=250)
        plt.close()

        plt.figure(figsize=(12,12)) 
        covered_both = dat[(dat.rna_umis >= 500) & (dat.dna_umis >=500)].groupby('sample_id')['rna_umis'].count().reset_index()
        covered_both.columns = ['sample', 'n_cells']
        covered_both.loc[:, 'assay'] = 'dna+rna'
        covered_rna = dat[(dat.rna_umis >=500)].groupby('sample_id')['rna_umis'].count().reset_index()
        covered_rna.columns = ['sample', 'n_cells']
        covered_rna.loc[:, 'assay'] = 'rna'
        covered_dna = dat[(dat.dna_umis >= 500)].groupby('sample_id')['dna_umis'].count().reset_index()
        covered_dna.columns = ['sample', 'n_cells']
        covered_dna.loc[:, 'assay'] = 'dna'
        
        ddf = pd.concat([covered_rna, covered_dna, covered_both])
        plt.figure(figsize=(12,8))
        sbn.barplot(x='sample', y='n_cells', hue='assay', data=ddf)
        plt.xticks(rotation=90);plt.ylabel('Number of covered cells');plt.xlabel('Sample')
        plt.tight_layout();plt.gca().set_rasterized(True)
        #plt.savefig(pdf, format='pdf', dpi=250)
        plt.close()

        
                


if __name__ == '__main__':
    main(get_args())
                


if __name__ == '__main__':
    main(get_args())
