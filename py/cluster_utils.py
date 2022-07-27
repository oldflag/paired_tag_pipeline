### Utilities for clustering-based QC

import numpy as np
import scanpy as scp
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import binned_statistic
import scipy as sp
from matplotlib import pyplot as plt
import seaborn as sbn
import warnings

PDF = None

MITO_GENES = ["ENSG00000210049", "ENSG00000211459", "ENSG00000210077", "ENSG00000210082", "ENSG00000209082",
              "ENSG00000198888", "ENSG00000210100", "ENSG00000210107", "ENSG00000210112", "ENSG00000198763",
              "ENSG00000210117", "ENSG00000210127", "ENSG00000210135", "ENSG00000210140", "ENSG00000210144",
              "ENSG00000198804", "ENSG00000210151", "ENSG00000210154", "ENSG00000198712", "ENSG00000210156",
              "ENSG00000228253", "ENSG00000198899", "ENSG00000198938", "ENSG00000210164", "ENSG00000198840",
              "ENSG00000210174", "ENSG00000212907", "ENSG00000198886", "ENSG00000210176", "ENSG00000210184",
              "ENSG00000210191", "ENSG00000198786", "ENSG00000198695", "ENSG00000210194", "ENSG00000198727",
              "ENSG00000210195", "ENSG00000210196", "MT-TF", "MT-RNR1", "MT-TV", "MT-RNR2", "MT-TL1", "MT-ND1",
              "MT-TI", "MT-TQ", "MT-TM", "MT-ND2", "MT-TW", "MT-TA", "MT-TN", "MT-TC", "MT-TY", "MT-CO1", "MT-TS1",
              "MT-TD", "MT-CO2", "MT-TK", "MT-ATP8", "MT-ATP6", "MT-CO3", "MT-TG", "MT-ND3", "MT-TR", "MT-ND4L",
              "MT-ND4", "MT-TH", "MT-TS2", "MT-TL2", "MT-ND5", "MT-ND6", "MT-TE", "MT-CYB", "MT-TT", "MT-TP",
              "ENSMUSG00000064336", "ENSMUSG00000064337", "ENSMUSG00000064338", "ENSMUSG00000064339",
              "ENSMUSG00000064340", "ENSMUSG00000064341", "ENSMUSG00000064342", "ENSMUSG00000064343",
              "ENSMUSG00000064344", "ENSMUSG00000064345", "ENSMUSG00000064346", "ENSMUSG00000064347",
              "ENSMUSG00000064348", "ENSMUSG00000064349", "ENSMUSG00000064350", "ENSMUSG00000064351",
              "ENSMUSG00000064352", "ENSMUSG00000064353", "ENSMUSG00000064354", "ENSMUSG00000064355",
              "ENSMUSG00000064356", "ENSMUSG00000064357", "ENSMUSG00000064358", "ENSMUSG00000064359",
              "ENSMUSG00000064360", "ENSMUSG00000064361", "ENSMUSG00000065947", "ENSMUSG00000064363",
              "ENSMUSG00000064364", "ENSMUSG00000064365", "ENSMUSG00000064366", "ENSMUSG00000064367",
              "ENSMUSG00000064368", "ENSMUSG00000064369", "ENSMUSG00000064370", "ENSMUSG00000064371",
              "ENSMUSG00000064372", "mt-Tf", "mt-Rnr1", "mt-Tv", "mt-Rnr2", "mt-Tl1", "mt-Nd1", "mt-Ti", "mt-Tq",
              "mt-Tm", "mt-Nd2", "mt-Tw", "mt-Ta", "mt-Tn", "mt-Tc", "mt-Ty", "mt-Co1", "mt-Ts1", "mt-Td", "mt-Co2",
              "mt-Tk", "mt-Atp8", "mt-Atp6", "mt-Co3", "mt-Tg", "mt-Nd3", "mt-Tr", "mt-Nd4l", "mt-Nd4", "mt-Th",
              "mt-Ts2", "mt-Tl2", "mt-Nd5", "mt-Nd6", "mt-Te", "mt-Cytb", "mt-Tt", "mt-Tp"]




def binary_jaccard(adata):
    """
    Binarize and compute the Jaccard matrix of an adata object

    Parameters
    ----------
    adata : An anndata object with an .X array

    Returns
    ----------
    Tuple of: Jaccard matrix, Expected matrix, Coverage vector
    """
    A = 1 * (adata.X > 0)
    and_ = A.dot(A.T)
    atot = and_.diagonal()
    or_ = atot[:, None] + atot - and_
    J = and_ / or_
    C = np.mean(A, axis=1).A1
    E = np.outer(C, C)/((C[np.newaxis,:] + C[:, np.newaxis]) - np.outer(C, C))
    return J.A, E, C


def lower_tri_sub(k, ix):
    """
    Return a subset (given by `ix`) of lower-triangular indices
    for a matrix of size `k x k`

    Parameters
    ----------
    k : the size of the matrix
    ix : An array of indexes (from 0 to k*(k-1)/2-1)

    Returns
    -------
    The tril_indices of k, subset to those in ix
    """
    tril_ix = np.tril_indices(k, -1)
    return tril_ix[0][ix], tril_ix[1][ix]


def normalize_jaccard(Jm, Em, nadj=500, winsor=99.0):
    """
    Normalize a Jaccard matrix given the expected overlap matrix

    Parameters
    ----------
    Jm : The Jaccard matrix
    Em : The expected overlap matrix (based on coverage)
    nadj : The number of datapoints to use for model fitting
    winsor : The winsorization level to apply to the normalized Jaccard matrix

    Returns
    -------
    A normalized Jaccard matrix

    """
    k = Jm.shape[0]
    ix = np.random.choice(int(k*(k-1)/2), nadj, replace=False)
    # fit the polynomial model
    poly = np.polyfit(
        Em[lower_tri_sub(k, ix)],
        Jm[lower_tri_sub(k, ix)],
        deg=2
    )
    # reset and apply
    Em[np.triu_indices(k, 0)] = 0
    Em[np.tril_indices(k, -1)] = np.polyval(poly,
                                            Em[np.tril_indices(k, -1)])
    Em += Em.T
    Em[np.diag_indices(k)] = 1
    J = Jm/Em

    # winsorize
    q = np.percentile(J[np.tril_indices(k, -1)], winsor)
    J[J > q] = q
    J[np.diag_indices(k)] = 1  # reset diagonal

    return J


def jaccard_diagnostics(dn_adata, n_approx=500, winsor=99.0, title=''):
    """\
    Compute the normalized Jaccard matrix for DNA data and produce
    the following diagnostic plots:

    [1]
     + Overall coverage histogram
     + Coverage by assay id
     + Un-normalized avg jaccard by coverage scatterplot
     + Normalized avg jaccard by coverage scatterplot
     + Normalized jaccard histogram
     + Jaccard diagonal boxplot

    [2]
     + Empaneled pairwise jaccard boxplot
    """
    J, E, C = binary_jaccard(dn_adata)
    Jn = normalize_jaccard(J, E, n_approx, winsor)
    D = Jn.sum(axis=1)
    assay_diagn = pd.DataFrame.from_dict({
        'assay_id': dn_adata.obs.assay_id.values,
        'bin_coverage': C,
        'jaccard_intensity': D
    })
    k = Jn.shape[0]
    ix = np.random.choice(int(k*(k-1)/2), 1000, replace=False)
    plt.figure()
    fig, axs = plt.subplots(2, 3, constrained_layout=True, figsize=(6*2.5, 4*2.5))
    axs[0, 0].hist(C, bins=50)
    axs[0, 0].set_xlabel('Cell mean coverage / bin')
    axs[0, 0].set_rasterized(True)
    sbn.boxplot(data=assay_diagn, x='assay_id', y='bin_coverage', ax=axs[0, 1])
    axs[0, 1].set_xticklabels(axs[0, 1].get_xticklabels(), rotation=90, fontsize=7)
    axs[0, 1].set_rasterized(True)
    axs[0, 2].scatter(C, (J.sum(axis=1)-1) / (2*(J.shape[0] - 1)), s=1)
    axs[0, 2].set_xlabel('Cell mean coverage / bin')
    axs[0, 2].set_ylabel('Cell average Jaccard')
    axs[0, 2].set_rasterized(True)
    axs[1, 0].scatter(C, (Jn.sum(axis=1)-1) / (2*(J.shape[0] - 1)), s=1)
    axs[1, 0].set_xlabel('Cell mean coverage / bin')
    axs[1, 0].set_ylabel('Cell average Corrected Jaccard')
    axs[1, 0].set_rasterized(True)
    axs[1, 1].hist(Jn[np.tril_indices(Jn.shape[0], -1)][ix])
    axs[1, 1].set_xlabel('Corrected Pairwise Jaccard')
    axs[1, 1].set_rasterized(True)
    sbn.boxplot(data=assay_diagn, x='assay_id', y='jaccard_intensity', ax=axs[1, 2])
    axs[1, 2].set_xticklabels(axs[1, 2].get_xticklabels(), rotation=90, fontsize=7)
    axs[1, 2].set_rasterized(True)

    plt.gca().set_rasterized(True)
    if PDF:
        plt.savefig(PDF, format='pdf', dpi=250);
        plt.close()

    # final Laplace-style normalization
    rs = 1 / np.sqrt(Jn.sum(axis=1))
    Jn = np.dot(np.diag(rs), np.dot(Jn, np.diag(rs)))

    ix2 = np.random.choice(int(k*(k-1)/2),
                           min(12000, int(k*(k-1)/2)))
    tril_ix = lower_tri_sub(k, ix2)
    bpdf = pd.DataFrame.from_dict({
        'row_index': tril_ix[0],
        'col_index': tril_ix[1],
        'assay1': dn_adata.obs.assay_id.values[tril_ix[0]],
        'assay2': dn_adata.obs.assay_id.values[tril_ix[1]],
        'jaccard': Jn[tril_ix]
    })
    plt.figure()

    assay_ids = sorted(dn_adata.obs.assay_id.unique().astype(str))
    if len(assay_ids) == 1:
        rr, cc = 1, 1
    elif len(assay_ids) == 2:
        rr, cc = 1, 2
    elif len(assay_ids) == 3:
        rr, cc, = 1, 3
    elif len(assay_ids) == 4:
        rr, cc = 2, 2
    else:
        rr, cc = 1 + int(np.sqrt(len(assay_ids))), 1+int(np.sqrt(len(assay_ids)))
    fig, axs = plt.subplots(rr, cc, constrained_layout=True, figsize=(6*2.5, 4*2.5))
    for i, aid in enumerate(assay_ids):
        c = i % cc
        r = min(int(i/cc), rr) - 1
        ddf = bpdf[bpdf.assay1 == aid]
        if len(assay_ids) < 4:
            sbn.boxplot(data=ddf, x='assay2', y='jaccard', ax=axs[c])
            axs[c].set_title(aid, size=8)
            axs[c].set_xticklabels(axs[c].get_xticklabels(), rotation=90, fontsize=7)
            axs[c].set_ylabel('Normalized Jaccard')
            axs[c].set_xlabel('Assay')
            axs[c].set_rasterized(True)
        else:
            sbn.boxplot(data=ddf, x='assay2', y='jaccard', ax=axs[r, c])
            axs[r, c].set_title(aid, size=8)
            axs[r, c].set_xticklabels(axs[r, c].get_xticklabels(), rotation=90, fontsize=7)
            axs[r, c].set_xlabel('Assay')
            axs[r, c].set_ylabel('Normalized Jaccard')
            axs[r, c].set_rasterized(True)
    plt.gca().set_rasterized(True)
    fig.suptitle(title)
    if PDF:
        plt.savefig(PDF, format='pdf', dpi=250)
        plt.close()

    return Jn


def do_dna_cluster_(dat_dna_analysis, npc, drop_first=True, resolution=0.5,
                    plot_title='DNA', n_neighbors=50, correct=None):
    """
    Perform DNA clustering (PC extraction, neighbor construction, UMAP, leiden)
    on a given jaccard-equipped adata object, and create the following
    empaneled diagnostic plot:
      + PCA 'elbow' plot
      + PC1/PC2 scatterplot, post-normalization, keyed by 'correct' or 'assay_id'
      + PC3/PC4 scatterplot, post-normalization, keyed by 'correct' or 'assay_id'
      + UMAP scatterplot, colored by leiden clusters
      + UMAP scatterplot, colored by 'correct' or 'assay_id'
      + t-SNE scatterplot, colored by leiden clusters



    Parameters
    ----------
    dat_dna_analysis : adata with a `.obsp['jaccard']` similarity matrix
    npc : the number of principal components to use
    drop_first : whether or not to drop the first principal component
    resolution : the resolution to use for Leiden clustering
    plot_title : the title to give the plots
    n_neighbors : the number of neighbors to use for UMAP
    correct : batch factor to correct for (run harmony); None for no correction

    Returns
    -------
    The input adata, updated with pca, umap, tsne slots; plus leiden clusters
    """
    # note that the `sparse` version is ordered small -> large
    U, s, Vt = sp.sparse.linalg.svds(dat_dna_analysis.obsp['jaccard'], k=npc + drop_first)
    if drop_first:
        U, Vt, s = U[:, -(npc+1):-1][:, ::-1], Vt[-(npc+1):-1, :][::-1], s[-(npc+1):-1][::-1]
    else:
        U, Vt, s = U[:, -npc:][:, ::-1], Vt[-npc, :][::-1, :], s[-npc:][::-1]

    dat_dna_analysis.obsm['X_pca'] = U

    if correct is not None:
        scp.external.pp.harmony_integrate(dat_dna_analysis, key=correct, basis='X_pca', adjusted_basis='X_pca_harmony',
                                          max_iter_harmony=150)
        rep_use='X_pca_harmony'
    else:
        rep_use='X_pca'

    dat_dna_analysis = scp.pp.neighbors(dat_dna_analysis, use_rep=rep_use, n_pcs=npc,
                                        copy=True,
                                        n_neighbors=n_neighbors)
    dat_dna_analysis = scp.tl.umap(dat_dna_analysis, copy=True)
    dat_dna_analysis = scp.tl.tsne(dat_dna_analysis, use_rep=rep_use, n_pcs=npc, copy=True)
    dat_dna_analysis = scp.tl.leiden(dat_dna_analysis, copy=True, resolution=resolution)

    fig, axs = plt.subplots(2, 3, constrained_layout=True, figsize=(6*2.5, 4*2.5))
    axs[0, 0].scatter(np.arange(s.shape[0]), s, color='black', s=2)
    axs[0, 0].set_xlabel('Eigenvalue #')
    axs[0, 0].set_ylabel('Eigenvalue')

    color_key = correct if correct is not None else 'assay_id'
    pca_data = pd.DataFrame.from_dict({
        'PC1': dat_dna_analysis.obsm[rep_use][:, 0],
        'PC2': dat_dna_analysis.obsm[rep_use][:, 1],
        'PC3': dat_dna_analysis.obsm[rep_use][:, 2],
        'PC4': dat_dna_analysis.obsm[rep_use][:, 3],
        'UMAP1': dat_dna_analysis.obsm['X_umap'][:, 0],
        'UMAP2': dat_dna_analysis.obsm['X_umap'][:, 1],
        'TSNE1': dat_dna_analysis.obsm['X_tsne'][:, 0],
        'TSNE2': dat_dna_analysis.obsm['X_tsne'][:, 1],
        color_key: dat_dna_analysis.obs.loc[:, color_key],
        'leiden': dat_dna_analysis.obs.leiden
    })
    sbn.scatterplot(data=pca_data, x='PC1', y='PC2', hue=color_key, legend=False, ax=axs[0, 1])
    sbn.scatterplot(data=pca_data, x='PC3', y='PC4', hue=color_key, legend=False, ax=axs[0, 2])
    sbn.scatterplot(data=pca_data, x='UMAP1', y='UMAP2', hue='leiden', legend=False, ax=axs[1, 0])
    sbn.scatterplot(data=pca_data, x='UMAP1', y='UMAP2', hue=color_key, legend=False, ax=axs[1, 1])
    sbn.scatterplot(data=pca_data, x='TSNE1', y='TSNE2', hue='leiden', legend=False, ax=axs[1, 2])

    for r in range(2):
        for c in range(3):
            axs[r, c].set_rasterized(True)

    plt.gca().set_rasterized(True)
    fig.suptitle(plot_title)
    if PDF:
        plt.savefig(PDF, format='pdf', dpi=250)
        plt.close()

    return dat_dna_analysis



def cluster_antibody_dna(dat_dna, antibody, n_pcs=15, drop_first=1, min_bins=300, min_umi=400, max_umi=3000,
                         min_cells=500, min_cells_bin=1, max_bin_pct=99.8, harmonize=None,
                         n_neighbors=50, resolution=0.5):
    """
    Subset a dna anndata object to a particular antibody and, if there are sufficiently many cells,
    perform clustering. Produce the following empaneled diagnostic plot:
       + Distribution of umi/bin; colored by autosome + x vs y + other
       + Filtered distribution of umi/bin
       + Count of # of cells / assay id
       + Count of # of cells / library id
       + UMI / cell distribution, pre-post filtering, by assay id
       + UMI / cell distribution, pre-post filtering, by library id

    """
    dat_dna = dat_dna[dat_dna.obs.antibody_name == antibody,:].copy()
    if dat_dna.shape[0] < min_cells:
        dat_dna.obs.loc[:,'leiden'] = np.nan
        return dat_dna
    dat_dna.raw = dat_dna
    dat_dna = dat_dna[(dat_dna.obs.feature_count >= min_bins) &
                      (dat_dna.obs.molecule_count >= min_umi) &
                      (dat_dna.obs.molecule_count <= max_umi)]

    fig, axs = plt.subplots(2, 3, constrained_layout=True, figsize=(6*2.5, 4*2.5))

    # now filter the features
    dat_dna.var.loc[:, 'counts'] = (1 * (dat_dna.X > 0)).sum(axis=0).A1
    dat_dna.var.loc[:, 'percentile'] = np.argsort(np.argsort(dat_dna.var.counts)) / dat_dna.var.shape[0]
    # autosome + X
    dat_dna.var.loc[:, 'contig'] = dat_dna.var.feature_name.map(lambda x: x.split('_')[0])
    want_contigs = {'chr%d' % s for s in range(1, 23)} | {'chrX'}
    dat_dna.var.loc[:, 'analysis'] = dat_dna.var.contig.isin(want_contigs).map(lambda s: 'retain' if s else 'remove')
    axs[0, 0].scatter(dat_dna.var[dat_dna.var.analysis == 'retain'].percentile.values,
                      dat_dna.var[dat_dna.var.analysis == 'retain'].counts.values - 0.25,  # slight offset
                      marker='.', label='Autosome + X')
    axs[0, 0].scatter(dat_dna.var[dat_dna.var.analysis == 'remove'].percentile.values,
                      dat_dna.var[dat_dna.var.analysis == 'remove'].counts.values + 0.25, # slight offset
                      marker='.', label='chrY + others')
    axs[0, 0].set_xlabel('Genomic bin (ranked by coverage)')
    axs[0, 0].set_ylabel('Bin coverage (number of cells with >=1 UMI)')
    axs[0, 0].set_rasterized(True)

    # subset bins
    dat_dna_analysis = dat_dna[:, (dat_dna.var.analysis == 'retain') &
                                  (dat_dna.var.percentile <= max_bin_pct) &
                                  (dat_dna.var.counts >= min_cells_bin)]


    axs[0, 1].scatter(dat_dna_analysis.var[dat_dna_analysis.var.analysis == 'retain'].percentile.values,
                      dat_dna_analysis.var[dat_dna_analysis.var.analysis == 'retain'].counts.values,
                      marker='.', label='Autosome + X')
    axs[0, 1].set_xlabel('Genomic bin (ranked by coverage)')
    axs[0, 1].set_ylabel('Bin coverage (number of cells in which bin is active)')
    axs[0, 1].set_rasterized(True)

    # recompute dna umi / cell
    dat_dna_analysis.obs.loc[:, 'filt_dna_umi'] = dat_dna_analysis.X.sum(axis=1).A1

    # num cells by assay / library
    cell_counts_assay = dat_dna_analysis.obs.groupby('assay_id')[['cell_id']].count().reset_index()
    sbn.barplot(data=cell_counts_assay, x='assay_id', y='cell_id', ax=axs[0, 2])
    axs[0, 2].set_xticklabels(axs[0, 2].get_xticklabels(), rotation=90, fontsize=7)
    axs[0, 2].set_xlabel('Assay')
    axs[0, 2].set_ylabel('# of Cells (pass-filter)')
    axs[0, 2].set_rasterized(True)

    cell_counts_library = dat_dna_analysis.obs.groupby('library_id')[['cell_id']].count().reset_index()
    sbn.barplot(data=cell_counts_library, x='library_id', y='cell_id', ax=axs[1,0])
    axs[1, 0].set_xticklabels(axs[1, 0].get_xticklabels(), rotation=90, fontsize=7)
    axs[1, 0].set_xlabel('Library')
    axs[1, 0].set_ylabel('# of Cells (pass-filter)')
    axs[1, 0].set_rasterized(True)

    # num umi / cell by assay / library
    sbn.boxplot(data=dat_dna_analysis.obs, x='assay_id', y='filt_dna_umi', ax=axs[1,1])
    axs[1, 1].set_xticklabels(axs[1, 1].get_xticklabels(), rotation=90, fontsize=7)
    axs[1, 1].set_xlabel('Assay')
    axs[1, 1].set_ylabel('# of UMI (per cell, post-bin filter)')
    axs[1, 1].set_rasterized(True)

    sbn.boxplot(data=dat_dna_analysis.obs, x='library_id', y='filt_dna_umi', ax=axs[1,2])
    axs[1, 2].set_xticklabels(axs[1, 2].get_xticklabels(), rotation=90, fontsize=7)
    axs[1, 2].set_xlabel('Library')
    axs[1, 2].set_ylabel('# of UMI (per cell, post-bin filter)')
    axs[1, 2].set_rasterized(True)

    fig.suptitle(antibody)
    plt.gca().set_rasterized(True)
    if PDF:
        plt.savefig(PDF, format='pdf', dpi=250)
        plt.close()

    dat_dna.obsp['jaccard'] = jaccard_diagnostics(dat_dna, n_approx=min_cells*4, title=antibody)
    dat_dna = do_dna_cluster_(dat_dna, n_pcs, drop_first, plot_title=antibody,
                              n_neighbors=n_neighbors, correct=harmonize,
                              resolution=resolution)
    return dat_dna


def find_modes1d(xdata, res=5000, bw=0.225):
    """
    Find peaks and troughs in a 1-d distribution of points by fitting a
    kernel density to it; and identifying points where the derivative
    moves +/- (peak) or -/+ (trough)

    """
    kde = sp.stats.gaussian_kde(xdata, bw_method=bw)
    x = np.arange(xdata.min(), xdata.max(), (xdata.max() - xdata.min()) / res)
    y = kde(x)
    plt.plot(x, y)
    d = np.hstack(([0], y[1:] - y[:-1]))
    p2n = np.hstack(([False], (d[:-1] > 0) & (d[1:] < 0)))
    n2p = np.hstack(([False], (d[:-1] < 0) & (d[1:] > 0)))
    return x[np.where(p2n)], x[np.where(n2p)]


def select_cells_assay_(anndata, assay_id, fallback_dna_umi, fallback_rna_umi, min_th=0.75, max_th=1.25):
    """
    Automatically select cells based on the dna/rna read distribution by searching for
    modes in a (typically) bimodal distribution.
    """
    assay_ix = np.where(anndata.obs.assay_id == assay_id)[0]
    r = np.sqrt(np.log10(anndata.obs.rna_reads.values[assay_ix]) ** 2 + \
                np.log10(anndata.obs.dna_reads.values[assay_ix]) ** 2)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        theta = np.arctan(np.log10(anndata.obs.dna_reads.values[assay_ix]) / \
                          np.log10(anndata.obs.rna_reads.values[assay_ix])) / (np.pi / 4)

    ix1 = np.where((theta >= min_th) & (theta <= max_th))[0]
    peaks, troughs = find_modes1d(r[ix1])
    lower, upper = peaks[(peaks > 1) & (peaks < 3)].mean(), peaks[(peaks > 4) & (peaks < 7)].mean()
    keep_vec = anndata.obs.loc[:, 'keep_cell_'].values
    if np.isnan(upper):
        print('Upper peak is too low for %s; falling back to UMI' % assay_id)
        dna_u, rna_u = anndata.obs.dna_umis.values[assay_ix], anndata.obs.rna_umis.values[assay_ix]
        ix2 = np.where((dna_u >= fallback_dna_umi) & (rna_u >= fallback_rna_umi))
        keep_idx = assay_ix[ix2]
    else:
        cut = troughs[(troughs > lower) & (troughs < upper)][-1]
        if np.isnan(cut) or np.isnan(lower) or np.isnan(upper):
            print('Cell distribution monotonic or unimodal for %s; falling back to UMI' % assay_id)
            keep_vec = anndata.obs.loc[:, 'keep_cell_'].values
            dna_u, rna_u = anndata.obs.dna_umis.values[assay_ix], anndata.obs.rna_umis.values[assay_ix]
            ix2 = np.where((dna_u >= fallback_dna_umi) & (rna_u >= fallback_rna_umi))
            keep_idx = assay_ix[ix2]
        else:
            cut2 = troughs[troughs > upper]
            if cut2.shape[0] > 0:
                ix2 = np.where((r[ix1] > cut) & (r[ix1] < cut2[-1]))
            else:
                ix2 = np.where(r[ix1] > cut)
            keep_idx = assay_ix[ix1][ix2]

    keep_vec[keep_idx] = True
    anndata.obs.loc[:, 'keep_cell_'] = keep_vec

    return anndata


def select_cells(anndata, fallback_dna_umi=350, fallback_rna_umi=350):
    """
    Automatically select cells based on the dna/rna read distribution by searching for
    modes in a (typically) bimodal distribution.

    Parameters
    ----------
    anndata : an anndata object produced by the PairedTag pipeline, with both dna and rna read + umi information
    fallback_dna_umi : The umi threshold to use if a good cut cannot be found
    fallback_rna_umi : The umi threshold to use if a good cut cannot be found

    Returns
    -------
    The input anndata object with a new `keep_cell_` column
    """
    print('Selecting cells...')
    anndata.obs.loc[:, 'keep_cell_'] = False  # set all to false (to generate the column)
    for assay_id in anndata.obs.assay_id.unique():
        anndata = select_cells_assay_(anndata, assay_id, fallback_dna_umi, fallback_rna_umi)

    assay_ids = anndata.obs.assay_id.unique()
    n_ids = len(assay_ids)
    rr, cc = int(np.sqrt(n_ids)) + 1, int(np.sqrt(n_ids)) + 1
    fig, axs = plt.subplots(rr, cc, figsize=(6 * 3, 6 * 3))
    for i, aid in enumerate(assay_ids):
        r = min(int(i / cc), cc) - 1
        c = i % cc
        dsub = anndata[anndata.obs.assay_id == aid]
        axs[r, c].scatter(dsub.obs.rna_reads, dsub.obs.dna_reads, s=1)
        dsub2 = dsub[dsub.obs.keep_cell_]
        axs[r, c].scatter(dsub2.obs.rna_reads, dsub2.obs.dna_reads, s=1)
        axs[r, c].set_xscale('log')
        axs[r, c].set_yscale('log')
        axs[r, c].set_title(aid, fontsize=8)
        axs[r, c].set_rasterized(True)

    plt.gca().set_rasterized(True)
    if PDF:
        plt.savefig(PDF, format='pdf', dpi=250)
        plt.close()

    return anndata


def cluster_pairedtag_dna(dat_dna, lim_features=250, lim_molecule=320, max_molecule=3000, n_pcs=10, drop_first=True,
                          min_cells=500, min_cells_bin=1, max_bin_excl=99.8, harmonize=None,
                          n_neighbors=50, resolution=0.5):
    """
    Cluster the dna (antibody) component of a PairedTag run. This will perform
    independent bin-based clustering of each antibody separately, and return a
    dictionary of {antibody -> clustered_adata}

    Parameters
    ----------
    dat_dna : An anndata object for dna, with .X, library_id, sample_id, antibody_name, and cell_id components
    lim_features : Minimum # bins a cell must have to include
    lim_molecule : Minimum # of molecules (UMI) a cell must have to
    max_molecule : Exclude cells with UMI above this threshold
    n_pcs : Number of principal components to use for clustering
    drop_first : Whether to drop the first PC. if n_pc=10 and drop_first=T, components 2-11 will be used rather than
        1-10
    min_cells : Minimum number of cells to perform clustering (else assign NA to the louvain component and return)
    min_cells_bin : Minimum number of 'active' cells in a bin to include it as a feature
    max_bin_excl : The (quantile) above which to exclude bins as potentially harboring mapping artifacts
    harmonize : The 'batch' feature column to use for harmony; else None for no batch correction
    n_neighbors : The number of neighbors to use for UMAP / Leiden
    resolution : The resolution to use for Leiden clustering

    Returns
    -------
    A dictionary of {antibody_name -> adata}, where the adata objects contain X_pca, X_umap, X_tsne, jaccard
    entries; as well as leiden clusters. If the adata is too small (< min_cells) then the adata object is
    unchanged, except that the 'leiden' component is set to NA.

    """
    antibs = dat_dna.obs.antibody_name.unique()
    if len(antibs) == 1:
        return do_dna_cluster_ab(dat_dna, antibs[0], n_pcs, drop_first,
                                 lim_features, lim_molecule, max_molecule, min_cells, min_cells_bin,
                                 max_bin_excl, harmonize, n_neighbors, resolution)
    else:
        return {ab: cluster_antibody_dna(dat_dna, ab, n_pcs, drop_first,
                                         lim_features, lim_molecule, max_molecule, min_cells, min_cells_bin,
                                         max_bin_excl, harmonize, n_neighbors, resolution)
                for ab in antibs if ab not in {'NA', 'Unused'}}


def scale_(u):
    return (u - np.mean(u)) / np.std(u)


def cluster_pairedtag_rna(obj, min_umi=350, max_umi=3500, min_cells_gene=10, n_pcs=10, res=0.5,
                          min_cells=500, n_genes=1000, n_neighbors=50, harmonize=None):
    """
    Cluster the RNA part of a PairedTag experiment and produce the following plots

      1. Empaneled (2 x 3)
        + Cell count per assay
        + Cell count per library
        + Cell UMI per assay
        + Cell UMI per library
        + % mtDNA per assay
        + % mtDNA per library

    2. Empaneled (2 x 3)
        + Variable gene plot
        + PCA Elbow plot
        + PC1 / PC2  (possibly post-harmony) colored by `harmonize` or assay
        + PC3 / PC4  (possibly post-harmony) colored by `harmonize` or assay
        + UMAP colored by leiden
        + Cell count per cluster [per assay or `harmonize`]

    3. Empaneled (2 x 3)
        + UMAP colored by leiden
        + UMAP colored by `harmonize` or assay
        + UMAP colored by antibody_name
        + tSNE colored by leiden
        + UMI per cluster
        + % mtDNA per cluster

    Parameters
    ----------
    obj : PairedTag anndata object containing mRNA counts, and RNA/DNA total UMI and total read counts
    min_umi : The minimum UMI to use for defining cells
    max_umi : The maximum UMI to use for defining cells (above this: likely doublet)
    min_cells_gene : The minimum number of cells to retain a gene
    n_pcs : The number of principal components to use for deriving embeddings
    res : The resolution to use for leiden clustering
    min_cells : The minimum number of cells to use for clustering (otherwise skip this entire function)
    n_genes : The number of variable genes to select
    n_neighbors : The number of nearest neighbors to use for Leiden / UMAP
    harmonize : A batch factor to adjust using Harmony (None for no batches)

    Returns
    -------
    The input object with PCA, (optionally harmonized) PCA, UMAP, and tSNE embeddings;
    as well as leiden clusters and % mitochoncrial UMI annotations.

    """
    obj.raw = obj
    print('Filtering RNA...')
    obj = obj[(obj.obs.molecule_count >= min_umi) & (obj.obs.molecule_count <= max_umi), :]
    obj.obs.loc[:, 'mito_UMI'] = obj[:, obj.var.feature_name.isin(MITO_GENES)].X.sum(axis=1)
    obj.obs.loc[:, 'pct_mito'] = 100 * obj.obs.mito_UMI / obj.obs.rna_umis
    scp.pp.filter_genes(obj, min_cells=min_cells_gene)

    if obj.shape[0] < min_cells:
        obj.obs.loc[:, 'leiden'] = np.nan
        return obj

    fig, axs = plt.subplots(2, 3, constrained_layout=True, figsize=(6 * 2.5, 4 * 2.5))
    sbn.countplot(data=obj.obs, x='assay_id', ax=axs[0, 0])
    axs[0, 0].set_xticklabels(axs[0, 0].get_xticklabels(), rotation=90, fontsize=7)
    axs[0, 0].set_xlabel('Assay ID')
    axs[0, 0].set_ylabel('# of cells (pass-filter)')
    axs[0, 0].set_rasterized(True)

    sbn.countplot(data=obj.obs, x='library_id', ax=axs[0, 1])
    axs[0, 1].set_xticklabels(axs[0, 1].get_xticklabels(), rotation=90, fontsize=7)
    axs[0, 1].set_xlabel('Library ID')
    axs[0, 1].set_ylabel('# of cells (pass-filter)')
    axs[0, 1].set_rasterized(True)

    sbn.boxplot(data=obj.obs, x='assay_id', y='rna_umis', ax=axs[0, 2])
    axs[0, 2].set_xticklabels(axs[0, 2].get_xticklabels(), rotation=90, fontsize=7)
    axs[0, 2].set_xlabel('Assay ID')
    axs[0, 2].set_ylabel('# of UMI (RNA)')
    axs[0, 2].set_rasterized(True)

    sbn.boxplot(data=obj.obs, x='library_id', y='rna_umis', ax=axs[1, 0])
    axs[1, 0].set_xticklabels(axs[1, 0].get_xticklabels(), rotation=90, fontsize=7)
    axs[1, 0].set_xlabel('Library ID')
    axs[1, 0].set_ylabel('# of UMI (RNA)')
    axs[1, 0].set_rasterized(True)

    sbn.boxplot(data=obj.obs, x='assay_id', y='pct_mito', ax=axs[1, 1])
    axs[1, 1].set_xticklabels(axs[1, 1].get_xticklabels(), rotation=90, fontsize=7)
    axs[1, 1].set_xlabel('Assay ID')
    axs[1, 1].set_ylabel('% Mitochondrial RNA')
    axs[1, 1].set_rasterized(True)

    sbn.boxplot(data=obj.obs, x='library_id', y='pct_mito', ax=axs[1, 2])
    axs[1, 2].set_xticklabels(axs[1, 2].get_xticklabels(), rotation=90, fontsize=7)
    axs[1, 2].set_xlabel('Library ID')
    axs[1, 2].set_ylabel('% Mitochondrial RNA')
    axs[1, 2].set_rasterized(True)

    plt.gca().set_rasterized(True)
    if PDF:
        plt.savefig(PDF, format='pdf', dpi=250)
        plt.close()
    
    print('Variable genes + PCA...')


    scp.pp.normalize_total(obj)
    scp.pp.log1p(obj)
    scp.pp.highly_variable_genes(obj, n_top_genes=n_genes, flavor='cell_ranger')
    fig, axs = plt.subplots(2, 3, constrained_layout=True, figsize=(6 * 2.5, 4 * 2.5))
    g1 = obj.var.highly_variable
    # variable genes need to be plotted manually on the axis (this is stupid)
    for lab, col, msk in zip(['highly variable genes', 'other genes'],
                         ['black', 'grey'], [g1, ~g1]):
        axs[0, 0].scatter(obj.var.means[msk], obj.var.dispersions[msk],
                          label=lab, c=col, s=1)
    axs[0, 0].set_xlabel('Mean expression of genes')
    axs[0, 0].set_ylabel('Dispersion of genes)')
    y_min, y_max = np.min(obj.var.dispersions), np.max(obj.var.dispersions)
    y_min = 0.95 * y_min if y_min > 0 else 1.05 * y_min
    axs[0, 0].set_ylim((y_min, y_max))
    axs[0, 0].set_rasterized(True)
    
    scp.pp.scale(obj)
    scp.tl.pca(obj)
    
    axs[0, 1].scatter(np.arange(obj.uns['pca']['variance_ratio'].shape[0]),
                            obj.uns['pca']['variance_ratio'], s=5, color='k')
    axs[0, 1].set_xlabel('Principal Component')
    axs[0, 1].set_ylabel('Variance Ratio')
    axs[0, 1].set_yscale('log')
    axs[0, 1].vlines(n_pcs, obj.uns['pca']['variance_ratio'].min(),
                 obj.uns['pca']['variance_ratio'].max(), linestyle='dashed',
                 color='red')
    axs[0, 1].set_rasterized(True)

    if harmonize is not None:
        scp.external.pp.harmony_integrate(obj, key=harmonize, basis='X_pca', adjusted_basis='X_pca_harmony',
                                         max_iter_harmony=150)
        rep_use = 'X_pca_harmony'
    else:
        rep_use = 'X_pca'

    print('Clustering + Leiden...')

    scp.pp.neighbors(obj, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep=rep_use)
    scp.tl.umap(obj)
    scp.tl.leiden(obj, resolution=res)
    scp.tl.tsne(obj)

    color_key = harmonize if harmonize is not None else 'assay_id'
    pca_data = pd.DataFrame.from_dict({
        'PC1': obj.obsm[rep_use][:, 0],
        'PC2': obj.obsm[rep_use][:, 1],
        'PC3': obj.obsm[rep_use][:, 2],
        'PC4': obj.obsm[rep_use][:, 3],
        'UMAP1': obj.obsm['X_umap'][:, 0],
        'UMAP2': obj.obsm['X_umap'][:, 1],
        'TSNE1': obj.obsm['X_tsne'][:, 0],
        'TSNE2': obj.obsm['X_tsne'][:, 1],
        color_key: obj.obs.loc[:, color_key],
        'antibody_name': obj.obs.antibody_name,
        'leiden': obj.obs.leiden
    })
    sbn.scatterplot(data=pca_data, x='PC1', y='PC2', hue=color_key, legend=False, ax=axs[0, 2], s=3)
    axs[0, 2].set_rasterized(True)
    sbn.scatterplot(data=pca_data, x='PC3', y='PC4', hue=color_key, legend=False, ax=axs[1, 0], s=3)
    axs[1, 0].set_rasterized(True)
    sbn.scatterplot(data=pca_data, x='UMAP1', y='UMAP2', hue='leiden', legend=False, ax=axs[1, 1], s=3)
    axs[1, 1].set_rasterized(True)
    sbn.countplot(data=obj.obs, x='leiden', ax=axs[1, 2])
    axs[1, 2].set_xlabel('Cluster number (Leiden)')
    axs[1, 2].set_ylabel('# of cells')
    axs[1, 2].set_rasterized(True)

    plt.gca().set_rasterized(True)
    if PDF:
        plt.savefig(PDF, format='pdf', dpi=250)
        plt.close()

    fig, axs = plt.subplots(2, 3, constrained_layout=True, figsize=(6 * 2.5, 4 * 2.5))
    sbn.scatterplot(data=pca_data, x='UMAP1', y='UMAP2', hue='leiden', legend=False, ax=axs[0, 0], s=3)
    axs[0, 0].set_rasterized(True)
    batch_hue = harmonize if harmonize is not None else 'assay_id'
    sbn.scatterplot(data=pca_data, x='UMAP1', y='UMAP2', hue=batch_hue, legend=False, ax=axs[0, 1], s=3)
    axs[0, 1].set_rasterized(True)
    sbn.scatterplot(data=pca_data, x='UMAP1', y='UMAP2', hue='antibody_name', legend=False, ax=axs[0, 2], s=3)
    axs[0, 2].set_rasterized(True)
    sbn.scatterplot(data=pca_data, x='TSNE1', y='TSNE2', hue='leiden', legend=False, ax=axs[1, 0], s=3)
    axs[1, 0].set_rasterized(True)
    sbn.boxplot(data=obj.obs, x='leiden', y='rna_umis', ax=axs[1, 1])
    axs[1, 1].set_xlabel('Cluster (Leiden)')
    axs[1, 1].set_ylabel('# of UMI (RNA)')
    axs[1, 1].set_rasterized(True)
    sbn.boxplot(data=obj.obs, x='leiden', y='pct_mito', ax=axs[1, 2])
    axs[1, 2].set_xlabel('Cluster (Leiden)')
    axs[1, 2].set_ylabel('% Mitochondrial RNA')
    axs[1, 2].set_rasterized(True)

    plt.gca().set_rasterized(True)
    if PDF:
        plt.savefig(PDF, format='pdf', dpi=250)
        plt.close()

    return obj

